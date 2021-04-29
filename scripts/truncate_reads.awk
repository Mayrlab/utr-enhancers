#!/usr/bin/gawk -f
BEGIN {
    if (length(len) == 0) {
	print "truncate_reads.awk: Variable 'len' must be defined!" > "/dev/stderr";
	exit 1;
    }
    printf("[INFO] Truncating reads to length %d\n", len) > "/dev/stderr";
}
{
    if (NR%2 == 1) {
	print $1;
    } else {
	print substr($0, 1, len);
    }
    if (NR%4000000 == 0) {
	printf("[INFO] %d reads processed\n", NR/4) > "/dev/stderr";
    }
}
END {
    printf("[INFO] Processed %d reads total\n", NR/4) > "/dev/stderr"
}
