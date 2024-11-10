import json
import boto3
import base64
import datetime
import time
import logging

logging.basicConfig(level=logging.DEBUG)
boto3.set_stream_logger("boto3.resources", logging.DEBUG)

# Initialize clients
s3 = boto3.client("s3")
dynamodb = boto3.resource("dynamodb")
table = dynamodb.Table("Jobs")
ec2_client = boto3.client("ec2")
ssm_client = boto3.client("ssm")


def lambda_handler(event, context):
    try:
        instance_id = None
        # Parse the JSON payload
        job_name = event["jobName"]
        email = event["email"]
        file_key = event["fileKey"]
        bucket_name = event["bucketName"]

        output_files = [
            f"s3://tmrnd-prototype/{file_key}",
            f"s3://tmrnd-prototype/{job_name}_diffused.pdb",
            f"s3://tmrnd-prototype/{job_name}_diffused.zarrmd",
        ]

        print("Logging job start")
        # Store job info in DynamoDB with timestamp
        store_job_info(job_name, output_files)

        print("Starting instance...")

        ami_id = "ami-0ed54081bb082582b"  # Replace with your AMI ID
        instance_type = (
            "g4dn.xlarge"  # Replace with your desired instance type
        )
        instance_id = start_instance(
            ami_id, instance_type
        )  # `key_name` can be None if using SSM only

        wait_for_instance(instance_id)
        print("Instance ready!")

        run_ssm_command(
            instance_id,
            f"/bin/bash -c 'aws s3 cp s3://tmrnd-prototype/{file_key} /root/inputs/{job_name}'",
        )

        print("Downloaded file to instance")
        main_command = "/bin/bash -c 'source ~/.bashrc && conda activate mda && python /root/analysis.py'"
        print("Running main command on instance")
        run_command_id = run_ssm_command(instance_id, main_command)

        print("Waiting for main command to complete...")
        while True:
            time.sleep(5)
            response = ssm_client.get_command_invocation(
                CommandId=run_command_id, InstanceId=instance_id
            )
            if response["Status"] == "Success":
                print("Command executed successfully. Output uploaded to S3.")
                break
            elif response["Status"] in ["Failed", "TimedOut", "Cancelled"]:
                print(f"Command failed with status: {response['Status']}")
                raise Exception(
                    f"SSM command execution failed with status: {response['Status']}"
                )

        print("Command complete, updating job status")

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
    finally:
        if instance_id is not None:
            terminate_instance(instance_id)


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


def start_instance(ami_id, instance_type):
    print("Attempting to call run_instances")
    try:
        response = ec2_client.run_instances(
            ImageId=ami_id,
            InstanceType=instance_type,
            MinCount=1,
            MaxCount=1,
            IamInstanceProfile={"Name": "S3_Admin"},
            # InstanceMarketOptions={"MarketType": "spot"},
            TagSpecifications=[
                {
                    "ResourceType": "instance",
                    "Tags": [{"Key": "Purpose", "Value": "SSM-Command-Test"}],
                }
            ],
        )
        print(response)
        instance_id = response["Instances"][0]["InstanceId"]
        print("Instance started")
        return instance_id
    except Exception as e:
        print(f"Failed to start EC2 instance: {e}")
        raise


def wait_for_instance(instance_id):
    # Wait until the instance is in 'running' state
    print(f"Waiting for instance {instance_id} to reach 'running' state...")
    ec2_client.get_waiter("instance_running").wait(InstanceIds=[instance_id])
    print("Instance is running. Waiting for SSM to be ready...")

    time.sleep(120)


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
