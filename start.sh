#/bin/bash/
set -e

trap 'echo "Error occurred at line $LINENO"' ERR
cmake -B build
make -C build

>./logfile.txt
if compgen -G "./tmp_file/*.fasta*" > /dev/null; then
    rm ./tmp_file/*.fasta*
fi
    # // data/truncated_sequences_1500.fasta
    # // data/Ecoli.fasta
    # //steps-512_0_NOT_U.fasta
    # //23s_rRNA.fasta
    # // chr1.fasta
    # data/Ecoli.fasta
time ./build/Crossfire -i ./data/mt1x.fasta -o ./out_res/test.out.fasta > test_log.txt

python3 ./out_res/SPscore.py --input ./out_res/test.out.fasta --match 0 --mismatch 1 --gap1 1 --gap2 1

