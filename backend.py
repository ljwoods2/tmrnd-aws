import json
import boto3
import base64

lambda_client = boto3.client("lambda")
s3 = boto3.client("s3")


def lambda_handler(event, context):
    try:
        # Parse the incoming request body
        body = json.loads(event["body"])
        job_name = body["jobName"]
        file_content_base64 = body["file"]

        # Decode the Base64 encoded file content
        file_content = base64.b64decode(file_content_base64)

        # Define S3 bucket and file key
        file_key = f"{job_name}_original.pdb"
        bucket_name = "tmrnd-prototype"

        # Upload file content to S3
        s3.put_object(
            Bucket=bucket_name,
            Key=file_key,
            Body=file_content,
            ContentType="application/octet-stream",
        )

        # Prepare payload for the second Lambda with S3 file path
        payload = {
            "jobName": job_name,
            "email": body["email"],
            "fileKey": file_key,
            "bucketName": bucket_name,
        }

        # Invoke the main Lambda asynchronously
        response = lambda_client.invoke(
            FunctionName="run_ec2",  # Replace with the ARN or name of the main Lambda
            InvocationType="Event",  # Asynchronous invocation
            Payload=json.dumps(payload),
        )

        # Immediately return a success response to the caller
        return {
            "statusCode": 200,
            "body": json.dumps(
                {
                    "message": "Job submitted successfully. Processing in background.",
                    "status": "Submitted",
                }
            ),
        }

    except Exception as e:
        # Handle any invocation errors
        return {
            "statusCode": 500,
            "body": json.dumps(
                {"error": str(e), "message": "Failed to submit job"}
            ),
        }
