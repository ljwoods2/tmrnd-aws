import json
import boto3
import base64
import datetime

s3 = boto3.client("s3")
dynamodb = boto3.resource("dynamodb")
table = dynamodb.Table("Jobs")


def lambda_handler(event, context):
    try:
        # Scan the entire table to get all jobs
        response = table.scan()
        items = response.get("Items", [])

        # Sort items by Timestamp from newest to oldest
        sorted_items = sorted(
            items, key=lambda x: x["Timestamp"], reverse=True
        )

        # Return a success response with the sorted items
        return {"statusCode": 200, "body": json.dumps(sorted_items)}

    except Exception as e:
        # Return an error
        return {
            "statusCode": 500,
            "body": json.dumps(
                {"error": str(e), "message": "Failed to query jobs"}
            ),
        }
