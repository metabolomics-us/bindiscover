from sqlalchemy import create_engine
from flask import Flask,request
from flask_restful import Api, Resource, reqparse
import json
from groupqueryresult import GroupQueryResult
import pandas as pd
import time

my_server='binvestigate-parker.czbqhgrlaqbf.us-west-2.rds.amazonaws.com'
my_database='postgres'
my_dialect='postgresql'
my_driver='psycopg2'
my_username='postgres'
my_password='S7kB93DIT46$'
my_port='5432'
my_connection=f'{my_dialect}+{my_driver}://{my_username}:{my_password}@{my_server}/{my_database}'
my_engine=create_engine(my_connection)#,echo=True)

class HGDAResource(Resource):

    def post(self):

        start=time.time()
        connection=my_engine.connect()
        end=time.time()

        multiple_metadata_panda=pd.DataFrame.from_dict(request.json['metadata_datatable'])
        my_GroupQueryResult=GroupQueryResult(multiple_metadata_panda,request.json['bin_type'],connection)

        my_engine.dispose()

        return json.dumps(my_GroupQueryResult.result_panda.to_dict('records'))