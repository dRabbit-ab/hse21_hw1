# hse21_hw1

## **Команды:**

Сначала создаю необходимые для работы папки
1. mkdir data
2. mkdir fastq
3. mkdir seqtk
4. mkdir trimmed
5. mkdir report
6. mkdir report/fastqc
7. mkdir report/multiqc
8. mkdir report/trimmed_fastqc
9. mkdir report/trimmed_fastqc

Далее создаем ссылки на необходимые файлы
11. cd fastq
12. ls -1 /usr/share/data-minor-bioinf/assembly/* | xargs -tI{} ln -s {}
13. cd ..

Выбираем случайные 5 миллионов чтений для PE и 1.5 для MP
15. seqtk sample -s119 fastq/oil_R1.fastq 5000000 > seqtk/pe1.fastq
16. seqtk sample -s119 fastq/oil_R2.fastq 5000000 > seqtk/pe2.fastq
17. seqtk sample -s119 fastq/oilMP_S4_L001_R1_001.fastq 1500000 > seqtk/mp1.fastq
18. seqtk sample -s119 fastq/oilMP_S4_L001_R2_001.fastq 1500000 > seqtk/mp2.fastq 

Оценим качество полученных чтений с помощью fastqc и multiqc 
19. cd seqtk
20. ls * | xargs -tI{} fastqc -o ~/data/report/fastqc {}
21. multiqc -o multiqc fastqc

С помощью программ platanus_trim и platanus_internal_trim подрежем чтения по качеству и удалим праймеры 
и переместим в нужную папку
22. platanus_trim seqtk/pe1.fastq seqtk/pe2.fastq
23. platanus_internal_trim seqtk/mp1.fastq seqtk/mp2.fastq
24. mv *trimmed ~/data/trimmed/

Оценим качество полученных чтений после обработки с помощью fastqc и multiqc
25. ls * | xargs -tI{} fastqc -o ~/data/report/trimmed_fastqc/ {}
26. multiqc -o trimmed_multiqc/ trimmed_fastqc/

Далее с помощью программы “platanus assemble” соберем контиги из подрезанных чтений
'''bash
mkdir contigs
'''



