BEGIN { FS=":" }
NF > 6 { print $7; } 
