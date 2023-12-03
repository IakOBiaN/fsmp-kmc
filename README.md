## Program FsMP_kMC

The program is designed for modeling adsorption layers using the FsMP_kMC method using numerical potentials.
You can read more about the method in the following works: [J. Phys. Chem. C 2021, 125, 50, 27853–27864](https://pubs.acs.org/doi/10.1021/acs.jpcc.1c09086?ref=pdf), [Phys. Chem. Chem. Phys., 2022,24, 26111-26123](https://pubs.rsc.org/en/content/articlelanding/2022/CP/D2CP03380A), and [Phys. Chem. Chem. Phys., 2023,25, 31352-31362](https://pubs.rsc.org/en/content/articlelanding/2023/CP/D3CP03955B).

## Instruction

To run the program, you first need to download or generate your potential and place it in the folder "potentials". You can download the available ones using the link: [DOWNLOAD NUMERICAL POTENTIALS](https://1drv.ms/f/s!AmyLqEdRe5EYgdkXdo7VUsFQxyMmng?e=6Vi3NS).
Then you need to take one of the existing ones or write your own configuration file and compile it. It is recommended to use the clang++ compiler. The program can then be launched to perform calculations.

```bash
clang++ -O3 terephthalic_acid_chain.cpp -o program.o
./program.o
```

## Программа FsMP_kMC
Программа предназначена для моделирования адсорбционнхы слоёв методом FsMP_kMC с использованием численных потенциалов.
Подробнее прочитать о методе можно в следующих работах: [J. Phys. Chem. C 2021, 125, 50, 27853–27864](https://pubs.acs.org/doi/10.1021/acs.jpcc.1c09086?ref=pdf), [Phys. Chem. Chem. Phys., 2022,24, 26111-26123](https://pubs.rsc.org/en/content/articlelanding/2022/CP/D2CP03380A) и [Phys. Chem. Chem. Phys., 2023,25, 31352-31362](https://pubs.rsc.org/en/content/articlelanding/2023/CP/D3CP03955B).

## Инструкция для запуска
Для запуска программы сначала необходимо скачать или сгенерировать свой потенциал и поместить его в папку "potentials". Скачать имеющиеся можно по ссылке: [СКАЧАТЬ ЧИСЛЕННЫЕ ПОТЕНЦИАЛЫ](https://1drv.ms/f/s!AmyLqEdRe5EYgdkXdo7VUsFQxyMmng?e=6Vi3NS).
Затем необходимо взять один из имеющихся или написать свой конфигурационный файл и провести его компиляцию. Рекомендуется использовать компилятор clang++. Затем программу можно запустить для проведения рассчётов.

```bash
clang++ -O3 terephthalic_acid_chain.cpp -o program.o
./program.o
```
