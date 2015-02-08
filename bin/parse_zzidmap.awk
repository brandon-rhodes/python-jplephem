BEGIN {
    print "target_names = ["
}

/^ *BLTCOD/ {
    n=$NF
}

/^ *BLTNAM/ {
    sub("[^']*", "")
    print "    (" n ", " $0 "),"
}

END {
    print "    ]"
}
