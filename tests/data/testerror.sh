#!/bin/bash

if [ "$1" == "1" ]; then
	echo "written to stdout"
fi

if [ "$2" == "1" ]; then
	echo "written to stderr" 1>&2;
fi

if [ "$3" == "1" ]; then
	exit 1
fi
