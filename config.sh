# tabix
export PATH=/home/yshira/bin/tabix-0.2.6:$PATH

# python
export PYTHONHOME=/usr/local/package/python2.7/2.7.8
export PYTHONPATH=~/local/lib/python2.7/site-packages
export PATH=${PYTHONHOME}/bin:${PATH}
# export LD_LIBRARY_PATH=${PYTHONHOME}/lib:${LD_LIBRARY_PATH}
export LD_LIBRARY_PATH=~/local/lib:${PYTHONHOME}/lib:${LD_LIBRARY_PATH}

check_error()
{

if [ $1 -ne 0 ]; then
    echo "FATAL ERROR: pipeline script"
    echo "ERROR CODE: $1"
    exit $1
fi
}


