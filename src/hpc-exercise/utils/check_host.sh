#/bin/sh

HNAME=`hostname`

if [ $HNAME = "cse" ]; then
	echo "!!!!!!!!!!!!!!!!!!!"
	echo "!!!!!!! cse !!!!!!!"
	echo "!!!!!!!!!!!!!!!!!!!"
else
	echo "local"
fi
