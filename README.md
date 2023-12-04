## Program FsMP_kMC

The program is designed for modeling adsorption layers using the FsMP_kMC method using numerical potentials.
You can read more about the method in the following works: [J. Phys. Chem. C 2021, 125, 50, 27853–27864](https://pubs.acs.org/doi/10.1021/acs.jpcc.1c09086?ref=pdf), [Phys. Chem. Chem. Phys., 2022,24, 26111-26123](https://pubs.rsc.org/en/content/articlelanding/2022/CP/D2CP03380A), and [Phys. Chem. Chem. Phys., 2023,25, 31352-31362](https://pubs.rsc.org/en/content/articlelanding/2023/CP/D3CP03955B).

## Instruction

* To run the program, you first need to download (or generate) potential and place it in the "potentials" folder. You can download the available numerical potentials using the link: [DOWNLOAD NUMERICAL POTENTIALS](https://1drv.ms/f/s!AmyLqEdRe5EYgdkXdo7VUsFQxyMmng?e=6Vi3NS).
* Then you need to go to the "configs" folder (cd configs).
* Take one of the available ones (for example, terephthalic_acid_chain.cpp) or write your own configuration file.
* Compile it. It is recommended to use the clang++ compiler (clang++ -O3 terephthalic_acid_chain.cpp -o program.o).
* The program can then be launched to perform calculations (./program.o).

```bash
cd configs
clang++ -O3 terephthalic_acid_chain.cpp -o program.o
./program.o
```

## Программа FsMP_kMC
Программа предназначена для моделирования адсорбционнхы слоёв методом FsMP_kMC с использованием численных потенциалов.
Подробнее прочитать о методе можно в следующих работах: [J. Phys. Chem. C 2021, 125, 50, 27853–27864](https://pubs.acs.org/doi/10.1021/acs.jpcc.1c09086?ref=pdf), [Phys. Chem. Chem. Phys., 2022,24, 26111-26123](https://pubs.rsc.org/en/content/articlelanding/2022/CP/D2CP03380A) и [Phys. Chem. Chem. Phys., 2023,25, 31352-31362](https://pubs.rsc.org/en/content/articlelanding/2023/CP/D3CP03955B).

## Инструкция для запуска
* Для запуска программы сначала необходимо скачать или сгенерировать свой потенциал и поместить его в папку "potentials". Скачать имеющиеся численные потенциалы можно по ссылке: [СКАЧАТЬ ЧИСЛЕННЫЕ ПОТЕНЦИАЛЫ](https://1drv.ms/f/s!AmyLqEdRe5EYgdkXdo7VUsFQxyMmng?e=6Vi3NS).
* Затем необходимо перейти в папку "configs" (cd configs).
* Взять один из имеющихся или написать свой конфигурационный файл. Например, terephthalic_acid_chain.cpp.
* Провести его компиляцию. Рекомендуется использовать компилятор clang++ (clang++ -O3 terephthalic_acid_chain.cpp -o program.o).
* Затем программу можно запустить для проведения рассчётов (./program.o).

```bash
cd configs
clang++ -O3 terephthalic_acid_chain.cpp -o program.o
./program.o
```
