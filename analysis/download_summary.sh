#!/usr/bin/env bash

if [ $# != 3 ]
then
    echo "Usage: $0 project instance_id task_name"
    exit 1
fi

DIR=$(cd `dirname $0`; pwd)
project=$1
instance=$2
taskname=$3

odpscmd --config=$DIR/../conf/odps_config.ini -e "http GET /projects/$project/instances/$instance?instancedetail&taskname=$taskname"
