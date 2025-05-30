#/bin/bash/

    # // data/truncated_sequences_1500.fasta
    # // data/Ecoli.fasta
    # //steps-512_0_NOT_U.fasta
    # //23s_rRNA.fasta
    # // chr1.fasta
    # data/Ecoli.fasta
    # Variola.fasta
time ./build/Crossfire -i ./data/steps-512_0_NOT_U.fasta -k 39 -c true -o ./out_res/test.out.fasta > test_log.txt


python3 ./out_res/SPscore.py --input ./out_res/test.out.fasta --match 0 --mismatch 1 --gap1 1 --gap2 1



