# 1. pattern {action}
# $0 any line with two digits
# $1 is the first field, $2 second field etc

#Wants to get out ID, date, time, open, lat, lon, pressure, cv, dp, rad, up, vp

BEGIN {counter = 0}

$0 ~ /Track/ {counter++;counter2=0}
$1 ~ /[0-9]\.[0-9]/ {
    counter2++
    printf("%6d %6d %08d %02d:00 %02d %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f\n", counter, counter2, $2,
$3/100, $6, $8, $9, $10, $11, $12, $13, $14, $15 )}

