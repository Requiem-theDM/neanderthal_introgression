#!/usr/bin/env python3

from sys import argv

lines = []

with open(f"{argv[1]}") as file:
    file.readline()
    for line in file:
        fields = line.split()
        lines.append(f"{fields[0]}\t{fields[1]}\t{fields[2]}\t{fields[3]}\t{fields[4]}\t{fields[5]}\t{fields[6]}\t{fields[7]}\t{fields[8]}\t{fields[9]}\t{argv[2]}")

for line in lines:
    print(line)