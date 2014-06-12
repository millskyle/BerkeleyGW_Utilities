while [ 1==1 ]; do 
   git add $@ 2>&1
   git commit -m "autocommit at `date '+%y-%m-%d-%H-%M-%S'`" > /dev/null 2>&1
   sleep 10
done
