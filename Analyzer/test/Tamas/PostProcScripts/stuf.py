import datetime
from sqlalchemy import TIMESTAMP

date_str = '2022-03-15'  # Your date string
date_obj = datetime.datetime.strptime(date_str, '%Y-%m-%d')  # Convert string to datetime object
timestamp_obj = TIMESTAMP(date_obj)  # Convert datetime object to TIMESTAMP object

