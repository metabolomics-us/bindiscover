from sqlalchemy import create_engine

from venntablequery import VennTableQuery
 
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


class VennTableResource(Resource):

    def post(self):
        '''
        '''

        print(request.json)

        dropdown_triplet_selection_value=request.json["dropdown_triplet_selection_value"]
        slider_percent_present_value=request.json["slider_percent_present_value"]
        toggle_average_true_value=bool(request.json["toggle_average_true_value"])
        radio_items_filter_value=request.json["radio_items_filter_value"]
        
        # print(toggle_average_true_value)
        # print(type(toggle_average_true_value))
        # print('------------------------')

        my_VennTableQuery=VennTableQuery()
        my_VennTableQuery.build_query(
            dropdown_triplet_selection_value,
            slider_percent_present_value,
            toggle_average_true_value,
            radio_items_filter_value
        )

        connection=my_engine.connect()
        temp_cursor=connection.execute(
            my_VennTableQuery.query
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
