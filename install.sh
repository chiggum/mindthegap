make -f  $(dirname $0)/makefile
SCRIPTPATH=$( cd $(dirname $0) ; pwd -P )
sudo ln -s  $SCRIPTPATH/bin/mindthegap /usr/bin/mindthegap