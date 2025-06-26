#!/bin/bash
# Ultra-simple coexistence detection

check_coexistence() {
    awk 'BEGIN{z=0;c=0} $2<0.001{z++;c=0} $2>=0.001{c++} END{
        if(z>5 && c>50) print "Coexistence"
        else if(z>5) print "Solid-only" 
        else print "Liquid-only"
    }' "${1:-density.dat}"
}

check_coexistence "$@"
