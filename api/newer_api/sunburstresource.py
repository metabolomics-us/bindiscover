from sqlalchemy import create_engine
 
from sunburstquery import SunburstQuery

from flask import Flask,request
from flask_restful import Api, Resource, reqparse
import json

my_server='binvestigate-parker.czbqhgrlaqbf.us-west-2.rds.amazonaws.com'
my_database='postgres'
my_dialect='postgresql'
my_driver='psycopg2'
my_username='postgres'
my_password='S7kB93DIT46$'
my_port='5432'
my_connection=f'{my_dialect}+{my_driver}://{my_username}:{my_password}@{my_server}/{my_database}'
my_engine=create_engine(my_connection)#,echo=True)

class SunburstResource(Resource):

    def post(self):
        '''
        '''


        compound=request.json['compound']


        my_SunburstQuery=SunburstQuery()
        my_SunburstQuery.build_query(
            compound
        )


        connection=my_engine.connect()
        temp_cursor=connection.execute(
            my_SunburstQuery.query
        )


        if (temp_cursor.rowcount <= 0):
            connection.close()
            #https://stackoverflow.com/questions/8645250/how-to-close-sqlalchemy-connection-in-mysql
            my_engine.dispose()
            print('row count of final result cursor less than 1')
            return 'fail'
        else:
            temp_result=json.dumps([dict(r) for r in temp_cursor])
            connection.close()
            #https://stackoverflow.com/questions/8645250/how-to-close-sqlalchemy-connection-in-mysql
            my_engine.dispose()
            #print(temp_result)
            return temp_result   
