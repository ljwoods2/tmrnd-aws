import json
import boto3
import base64
import datetime

# Initialize the S3 client
s3 = boto3.client('s3')

def lambda_handler(event, context):
    try:
        # Parse the JSON payload
        body = json.loads(event['body'])
        email = body['email']
        job_name = body['jobName']
        file_content_base64 = body['file']
        
        # Decode the Base64 encoded file content
        file_content = base64.b64decode(file_content_base64)
        
        # Define the S3 bucket and file key
        bucket_name = 'job-files-bucket'  # Replace with your S3 bucket name
        file_key = f"{job_name}-{datetime.datetime.now().isoformat()}.pdb"  # Unique file name with timestamp

        # Upload the file to S3
        s3.put_object(
            Bucket=bucket_name,
            Key=file_key,
            Body=file_content,
            ContentType='application/octet-stream'
        )

        # Return a success response
        return {
            'statusCode': 200,
            'body': json.dumps({
                'message': 'Job submitted successfully',
                'email': email,
                'jobName': job_name,
                'fileKey': file_key
            })
        }

    except Exception as e:
        # Return an error response
        return {
            'statusCode': 500,
            'body': json.dumps({
                'error': str(e),
                'message': 'Failed to submit job'
            })
        }
