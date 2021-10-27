# Домашнее задание 1

## Обязательная часть

Сначала создаю необходимые для работы папки
```bash
mkdir data
mkdir fastq
mkdir seqtk
mkdir trimmed
mkdir report
mkdir report/fastqc
mkdir report/multiqc
mkdir report/trimmed_fastqc
mkdir report/trimmed_fastqc
```

Далее создаем ссылки на необходимые файлы
```bash
cd fastq
ls -1 /usr/share/data-minor-bioinf/assembly/* | xargs -tI{} ln -s {}
cd ..
```

Выбираем случайные 5 миллионов чтений для PE и 1.5 для MP
```bash
seqtk sample -s119 fastq/oil_R1.fastq 5000000 > seqtk/pe1.fastq
seqtk sample -s119 fastq/oil_R2.fastq 5000000 > seqtk/pe2.fastq
seqtk sample -s119 fastq/oilMP_S4_L001_R1_001.fastq 1500000 > seqtk/mp1.fastq
seqtk sample -s119 fastq/oilMP_S4_L001_R2_001.fastq 1500000 > seqtk/mp2.fastq 
```

Оценим качество полученных чтений с помощью fastqc и multiqc
```bash
cd seqtk
ls * | xargs -tI{} fastqc -o ~/data/report/fastqc {}
multiqc -o multiqc fastqc
```

С помощью программ platanus_trim и platanus_internal_trim подрежем чтения по качеству и удалим праймеры 
и переместим в нужную папку
```bash
platanus_trim seqtk/pe1.fastq seqtk/pe2.fastq
platanus_internal_trim seqtk/mp1.fastq seqtk/mp2.fastq
mv *trimmed ~/data/trimmed/
```

Оценим качество полученных чтений после обработки с помощью fastqc и multiqc
```bash
ls * | xargs -tI{} fastqc -o ~/data/report/trimmed_fastqc/ {}
multiqc -o trimmed_multiqc/ trimmed_fastqc/
```

Далее с помощью программы “platanus assemble” соберем контиги из подрезанных чтений
```bash
mkdir contigs
platanus assemble -o Poil -t 2 -m 20 -f ~/data/trimmed/pe1.fastq.trimmed ~/data/trimmed/pe2.fastq.trimmed 2> assemble.log
```

Соберем скаффолды из полученных контигов
```bash
mkdir scaf
platanus scaffold -o Poil -t 2 -c ~/data/contigs/Poil_contig.fa -IP1 ~/data/trimmed/pe1.fastq.trimmed ~/data/trimmed/pe2.fastq.trimmed -OP2 ~/data/trimmed/mp1.fastq.int_trimmed ~/data/trimmed/mp2.fastq.int_trimmed 2> scaffold.log
```

Уменьшим количество гэпов с помощью утилиты platanus gap_close
```bash
mkdir gap_close
platanus gap_close -o Poil -t 2 -c ~/data/scaf/Poil_scaffold.fa -IP1 ~/data/trimmed/pe1.fastq.trimmed ~/data/trimmed/pe2.fastq.trimmed -OP2 ~/data/trimmed/mp1.fastq.int_trimmed ~/data/trimmed/mp2.fastq.int_trimmed 2> gapclose.log
```

## Необязательная часть
