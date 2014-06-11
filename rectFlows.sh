#!/bin/bash

nice -n7 ./a.out 50 1 50 &
nice -n7 ./a.out 50 1 70 &
nice -ny ./a.out 50 1 100 &
nice -n7 ./a.out 50 1 150 &
nice -n7 ./a.out 50 1 200 &

nice -n7 ./a.out 100 1 50 &
nice -n7 ./a.out 100 1 70 &
nice -ny ./a.out 100 1 100 &
nice -n7 ./a.out 100 1 150 &
nice -n7 ./a.out 100 1 200 &

nice -n7 ./a.out 200 1 50 &
nice -n7 ./a.out 200 1 70 &
nice -ny ./a.out 200 1 100 &
nice -n7 ./a.out 200 1 150 &
nice -n7 ./a.out 200 1 200 &
