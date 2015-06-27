# tabix
export PATH=/home/yshira/bin/tabix-0.2.6:$PATH

# python
export PYTHONHOME=/usr/local/package/python2.7/2.7.8
export PYTHONPATH=${PYTHONHOME}/lib/python2.7/site-packages:~/local/lib/python2.7/site-packages
export PATH=${PYTHONHOME}/bin:${PATH}
# export LD_LIBRARY_PATH=${PYTHONHOME}/lib:${LD_LIBRARY_PATH}
export LD_LIBRARY_PATH=~/local/lib:${PYTHONHOME}/lib:${LD_LIBRARY_PATH}

# reference
REFERENCE=/home/ogawaprj/ngs/ref/GRCh37-lite_PCAWG_AB513134_bwa-0.7.10/GRCh37-lite_PCAWG_AB513134.fa

check_error()
{

if [ $1 -ne 0 ]; then
    echo "FATAL ERROR: pipeline script"
    echo "ERROR CODE: $1"
    exit $1
fi
}


