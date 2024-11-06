import json
import boto3
import base64
import datetime
import time

# Initialize clients
s3 = boto3.client("s3")
dynamodb = boto3.resource("dynamodb")
table = dynamodb.Table("Jobs")
ec2_client = boto3.client("ec2")
ssm_client = boto3.client("ssm")


def lambda_handler(event, context):
    try:
        # Parse the JSON payload
        body = json.loads(event["body"])
        email = body["email"]
        job_name = body["jobName"]
        file_content_base64 = body["file"]

        # Decode the Base64 encoded file content
        file_content = base64.b64decode(file_content_base64)

        # Define the S3 bucket and file key
        bucket_name = "tmrnd-prototype"  # Replace with your S3 bucket name
        file_key = f"{job_name}-{datetime.datetime.now().isoformat()}.pdb"  # Unique file name with timestamp
        output_files = [file_key]

        # Store job info in DynamoDB with timestamp
        store_job_info(job_name, output_files)

        # Upload the file to S3
        s3.put_object(
            Bucket=bucket_name,
            Key=file_key,
            Body=file_content,
            ContentType="application/octet-stream",
        )

        # Update job status to completed and set last updated timestamp
        update_job_status(job_name, "COMPLETED")

        # Return a success response
        return {
            "statusCode": 200,
            "body": json.dumps(
                {
                    "message": "Job submitted successfully",
                    "email": email,
                    "jobName": job_name,
                    "fileKey": file_key,
                }
            ),
        }

    except Exception as e:
        # Attempt to update status to FAILED if an error occurs
        try:
            update_job_status(job_name, "FAILED")
        except:
            pass
        finally:
            # Return an error response
            return {
                "statusCode": 500,
                "body": json.dumps(
                    {"error": str(e), "message": "Failed to submit job"}
                ),
            }


def store_job_info(job_id, output_files):
    """Store job information with initial status, output files, and timestamp in DynamoDB."""
    table.put_item(
        Item={
            "JobID": job_id,
            "Status": "IN_PROGRESS",
            "OutputFiles": output_files,  # List of S3 URLs for output files
            "Timestamp": datetime.datetime.now().isoformat(),  # Initial creation timestamp
        }
    )


def update_job_status(job_id, new_status):
    """Update the job status in DynamoDB."""
    response = table.update_item(
        Key={"JobID": job_id},
        UpdateExpression="SET #st = :status",
        ExpressionAttributeNames={"#st": "Status"},
        ExpressionAttributeValues={
            ":status": new_status,
        },
        ReturnValues="UPDATED_NEW",
    )
    return response


def start_instance(ami_id, instance_type, key_name):
    response = ec2_client.run_instances(
        ImageId=ami_id,
        InstanceType=instance_type,
        MinCount=1,
        MaxCount=1,
        KeyName=key_name,
        IamInstanceProfile={"Name": "tmrnd"},
        InstanceMarketOptions={"MarketType": "spot"},
        TagSpecifications=[
            {
                "ResourceType": "instance",
                "Tags": [{"Key": "Purpose", "Value": "SSM-Command-Test"}],
            }
        ],
    )
    instance_id = response["Instances"][0]["InstanceId"]
    return instance_id


def wait_for_instance(instance_id):
    waiter = ec2_client.get_waiter("instance_running")
    waiter.wait(InstanceIds=[instance_id])


def download_file_to_instance(instance_id, bucket_name, s3_key, local_path):
    """Send an SSM command to download a file from S3 onto the instance."""
    command = f"aws s3 cp s3://{bucket_name}/{s3_key} {local_path}"
    response = ssm_client.send_command(
        InstanceIds=[instance_id],
        DocumentName="AWS-RunShellScript",
        Parameters={"commands": [command]},
    )
    command_id = response["Command"]["CommandId"]
    return command_id


def run_ssm_command(instance_id, command):
    response = ssm_client.send_command(
        InstanceIds=[instance_id],
        DocumentName="AWS-RunShellScript",
        Parameters={"commands": [command]},
    )
    command_id = response["Command"]["CommandId"]
    return command_id


def upload_output_to_s3(instance_id, command_id, bucket_name, s3_key):
    print(f"Waiting for command {command_id} to complete...")
    while True:
        time.sleep(5)
        response = ssm_client.get_command_invocation(
            CommandId=command_id, InstanceId=instance_id
        )
        if response["Status"] == "Success":
            print("Command executed successfully. Uploading output to S3.")
            # Retrieve output from instance and upload to S3
            output = response["StandardOutputContent"]
            s3.put_object(Bucket=bucket_name, Key=s3_key, Body=output)
            print(f"Output uploaded to s3://{bucket_name}/{s3_key}")
            break
        elif response["Status"] in ["Failed", "TimedOut", "Cancelled"]:
            print(f"Command failed with status: {response['Status']}")
            raise Exception(
                f"SSM command execution failed with status: {response['Status']}"
            )


def terminate_instance(instance_id):
    print(f"Terminating instance {instance_id}...")
    ec2_client.terminate_instances(InstanceIds=[instance_id])
    print(f"Instance {instance_id} terminated.")
