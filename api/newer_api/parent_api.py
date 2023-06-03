from pprint import pprint
from flask import Flask,request
from flask_restful import Api, Resource, reqparse
import json

from sqlalchemy import create_engine
from sqlalchemy import Table, String
from sqlalchemy.dialects import postgresql

from venntableresource import VennTableResource
from sunburstresource import SunburstResource
from leafresource import LeafResource
from leafmetadataresource import LeafMetadataResource
from hgdametadataresource import HGDAMetadataResource
from another_name import HGDAResource
from binresource import BinResource
from treeresourceextra import TreeResource
from treemetadataresource import TreeMetadataResource
#there is something really really weird going on where it would nto import from treeresoure or hgdaresource so i renamed them both
#and it worked

app=Flask(__name__)
api=Api(app)

my_server='binvestigate-parker.czbqhgrlaqbf.us-west-2.rds.amazonaws.com'
my_database='postgres'
my_dialect='postgresql'
my_driver='psycopg2'
my_username='postgres'
my_password='S7kB93DIT46$'
my_port='5432'
my_connection=f'{my_dialect}+{my_driver}://{my_username}:{my_password}@{my_server}/{my_database}'
my_engine=create_engine(my_connection)#,echo=True)


api.add_resource(VennTableResource,'/venntableresource/')
api.add_resource(SunburstResource,'/sunburstresource/')
api.add_resource(LeafMetadataResource,'/leafmetadataresource/')
api.add_resource(HGDAMetadataResource,'/hgdametadataresource/')
api.add_resource(LeafResource,'/leafresource/')
api.add_resource(HGDAResource,'/hgdaresource/')
api.add_resource(BinResource,'/binresource/')
api.add_resource(TreeMetadataResource,'/treemetadataresource/')
api.add_resource(TreeResource,'/treeresource/')


if __name__ == '__main__':
    app.run(debug=False,port=4999,host='0.0.0.0')
    #app.run(debug=True,port=4999)
