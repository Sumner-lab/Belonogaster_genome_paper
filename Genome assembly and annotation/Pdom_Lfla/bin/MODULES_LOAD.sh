#!/bin/bash

to_load=$1
do_stack=$2

if [ "$do_stack" != "stack" ]
then
	if [ ! -z "$do_stack" ]
	then
      echo "2nd argument value must be \"\" or \"stack\"."
    fi
do_stack="no"
fi

module purge

module_file=$HOME/config_files/module_list.cfg

module_to_load=$(grep "^${to_load}:" $module_file | cut -f 2 -d ":" | sed "s/ //g")

. $module_to_load $do_stack

if [ $? -eq 0 ]
then
    echo "Module loaded:$module_to_load"
else
    echo "ERROR loading module ${to_load}: something failed"
	echo "Retrying debug"
	. $module_to_load $do_stack debug
	if [ $? -eq 0 ]
    then
      echo "Module loaded:$module_to_load"
	else
	  echo "Module ${to_load} cannot be loaded"
	  exit
	fi
fi

