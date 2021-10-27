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
mkdir report/trimmed_multiqc
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

![image](https://user-images.githubusercontent.com/79662580/139108683-59ffc23a-8fb2-4527-8185-5b3100a78b9f.png)
![fastqc_adapter_content_plot](https://user-images.githubusercontent.com/79662580/139109521-867bdb7a-e7ed-46c9-af2f-382729ef36f0.png)
![fastqc_per_base_sequence_quality_plot](https://user-images.githubusercontent.com/79662580/139109528-b290503e-3aeb-47c2-a660-a4ef7ddfc0f7.png)
![fastqc_per_sequence_quality_scores_plot](https://user-images.githubusercontent.com/79662580/139109535-8ba60fff-d285-484d-b29f-52e645d20b92.png)
![fastqc_sequence_counts_plot](https://user-images.githubusercontent.com/79662580/139109546-1997a18d-794f-4a05-9c8e-325b337fe310.png)
По статистике видно, что существуют некоторые проблемы связанные с наличием адапторов, что портит качество чтений

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

![image](https://user-images.githubusercontent.com/79662580/139109215-a30349c5-8161-414e-b5a5-d3bc70f345ba.png)
![fastqc_adapter_content_plot](https://user-images.githubusercontent.com/79662580/139116346-66622ab5-f6f3-4bf0-ae21-2a1c25ae223f.png)
![fastqc_per_base_sequence_quality_plot](https://user-images.githubusercontent.com/79662580/139116357-1e1e77d0-2601-461f-a93e-55d186a4edf4.png)
![fastqc_per_sequence_quality_scores_plot](https://user-images.githubusercontent.com/79662580/139116365-480aea2b-452b-4d35-82c7-23e9622ba307.png)
![fastqc_sequence_counts_plot](https://user-images.githubusercontent.com/79662580/139116373-3428420e-d529-4fb1-8486-8fd9ab4a36c7.png)
После обработки видно, что качество заметно улучшилось, но уменьшилась длина чтений, так как были удалены адаптеры

Далее с помощью программы platanus assemble соберем контиги из подрезанных чтений
```bash
mkdir contigs
platanus assemble -o Poil -t 2 -m 20 -f ~/data/trimmed/pe1.fastq.trimmed ~/data/trimmed/pe2.fastq.trimmed 2> assemble.log
```

Сразу напишем код, который нам понадобится для анализа (анализ гэпов, контигов, загрузка данных)
```python
def get_data(path):
    data = list()
    with open(path, 'r') as file:
        for line in file:
            if line[0] == '>':
                data.append('')
            else:
                data[-1] += line[:-1]
        return data
    
    
def n_50(data, length):
    data = sorted([len(item) for item in data], reverse=True)
    s = 0
    counter = 0
    while True:
        s += data[counter]
        if s / length >= 0.5:
            return data[counter]
        counter += 1
    

def get_stats(data):
    result = dict()
    result['Contig_number'] = len(data)
    result['General_length'] = sum([len(item) for item in data])
    result['Max_length'] = max([len(item) for item in data])
    result['N50'] = n_50(data, result['General_length'])
    return result


def gap_analysis(data):
    data = sorted(data, reverse=True, key=lambda x: len(x))[0]
    print('success_1')
    s = 0
    quantity = 0
    for i in range(len(data)):
        if data[i] == 'N':
            s += 1
            
    if data[0] == 'N':
        inside = True
        quantity += 1
    else:
        inside = False
        
    for i in range(len(data)):
        if data[i] == 'N':
            if not inside:
                quantity += 1
            inside = True
        else:
            inside = False
    print('success_3')
    return quantity, s
```

Посмотрим на результаты
![image](https://user-images.githubusercontent.com/79662580/139112049-d22e88b4-9179-49fb-9618-ece491e835cb.png)

Соберем скаффолды из полученных контигов с помощью утилиты platanus scaffold
```bash
mkdir scaf
platanus scaffold -o Poil -t 2 -c ~/data/contigs/Poil_contig.fa -IP1 ~/data/trimmed/pe1.fastq.trimmed ~/data/trimmed/pe2.fastq.trimmed -OP2 ~/data/trimmed/mp1.fastq.int_trimmed ~/data/trimmed/mp2.fastq.int_trimmed 2> scaffold.log
```

Теперь посмотрим, какая статистика получается для скаффолдов
![image](https://user-images.githubusercontent.com/79662580/139112291-06589d52-86c6-430c-8004-57dafdc0ce86.png)
Видим, что количество скаффолдов значительно меньше, чем количество контигов, а самый длинный имеет почти общую длинну последовательности, что означает, 
что это фактически и есть собранный геном

Теперь обратимся к гэпам
![image](https://user-images.githubusercontent.com/79662580/139112751-f4172002-2361-4e22-8a74-0a1e75a94779.png)

Уменьшим количество гэпов с помощью утилиты platanus gap_close
```bash
mkdir gap_close
platanus gap_close -o Poil -t 2 -c ~/data/scaf/Poil_scaffold.fa -IP1 ~/data/trimmed/pe1.fastq.trimmed ~/data/trimmed/pe2.fastq.trimmed -OP2 ~/data/trimmed/mp1.fastq.int_trimmed ~/data/trimmed/mp2.fastq.int_trimmed 2> gapclose.log
```
Посмотрим, что изменилось
![image](https://user-images.githubusercontent.com/79662580/139112944-e0d8eb45-ca10-40a2-bdd7-48af3ab2a48a.png)
Количество гэпов значительно уменшилось, как и их длинна

## Необязательная часть
Попробуем значительно уменьшить количество чтений, например, до 1 миллиона для PE и 200 тысяч для MP и посмотрим на результат
![image](https://user-images.githubusercontent.com/79662580/139117703-9fc8d8a3-1af9-4448-b869-8bf75efbd7c2.png)
![fastqc_adapter_content_plot](https://user-images.githubusercontent.com/79662580/139117774-b3c5e70a-e343-41c9-aaae-7e26b0b4dd38.png)
![fastqc_per_base_sequence_quality_plot](https://user-images.githubusercontent.com/79662580/139117791-e00a9306-3538-4f5f-890c-24e3657addae.png)
![fastqc_per_sequence_quality_scores_plot](https://user-images.githubusercontent.com/79662580/139117806-bb67525b-309e-4081-93f3-4b270d8b3268.png)
![fastqc_sequence_counts_plot](https://user-images.githubusercontent.com/79662580/139117822-0b44abf7-5ee8-42f8-96a6-27516767ca08.png)

Теперь обработаем адаптеры и посмотрим на статистику
![image](https://user-images.githubusercontent.com/79662580/139119296-dc16f839-c444-4389-bb35-18052a43d930.png)
![fastqc_adapter_content_plot](https://user-images.githubusercontent.com/79662580/139119407-beaff209-7d23-41f4-b66a-0b526455fa31.png)
![fastqc_per_base_sequence_quality_plot](https://user-images.githubusercontent.com/79662580/139119413-983759e3-1141-44af-afa9-715d293986b7.png)
![fastqc_per_sequence_quality_scores_plot](https://user-images.githubusercontent.com/79662580/139119444-ccacc861-d374-465c-ade1-5a61af3554ec.png)
![fastqc_sequence_counts_plot](https://user-images.githubusercontent.com/79662580/139119456-c43d6eee-789e-4e5b-9389-f095ed9a7526.png)
Удаление адаптеров все так же хорошо исправляет чтения, так же было и в первом примере, ничего особо не поменялось

Перейдем к сборке контигов и скаффолдов
![image](https://user-images.githubusercontent.com/79662580/139121655-d5214880-d6b6-4b72-9e54-53cea4f0dcca.png)
Скаффолдов все так же меньше, чем контигов, как и в первом примере, но вот их число выросло, а длинна самого большого сильно уменьшилась, это значительно снижает качество выполненой работы. Количество гэпов уменьшилось, но gap_close сработал менее эффективно, нежели в первом примере. Вывод - качество сборки значительно упало
