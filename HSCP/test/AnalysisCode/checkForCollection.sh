collection=$1
dir1="../"
dir2="../../plugins"
echo "#### List of calls for Collection $collection"
echo "======================================================================================"
find $dir1 -maxdepth 3 -type f -exec grep $2 --color $collection {} +
find $dir2 -maxdepth 3 -type f -exec grep $2 --color $collection {} +
echo "======================================================================================"
