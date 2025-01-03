#/bin/bash/
set -e

trap 'echo "Error occurred at line $LINENO"' ERR
cmake -B build
make -C build

>./logfile.txt
if compgen -G "./tmp_file/*.fasta*" > /dev/null; then
    rm ./tmp_file/*.fasta*
fi

./build/Crossfire -i ./data/test.fasta -o ./out_res/test.out.fasta

python3 ./out_res/SPscore.py --input ./out_res/test.out.fasta --match 0 --mismatch 1 --gap1 1 --gap2 1

