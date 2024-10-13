# Инструкция по установке SLEPc

### Оглавление
1. [Arch Linux](#arch-linux)
   - [Установка PETSc](#установка-petsc)
   - [Установка SLEPc](#установка-slepc)
   - [Проверка установки](#проверка-установки)
2. [Windows](#windows)
3. [Полезные ссылки](#полезные-ссылки)

## Arch Linux
1. Склонируйте репозитории SLEPc и PETSc с ветками release:
```
git clone -b release https://gitlab.com/petsc/petsc.git petsc
git clone -b release git@github.com:slepc/slepc.git
```

### Установка PETSc

2. Задайте переменные окружения:
```
export PETSC_DIR=/home/jess/develop/slepc/petsc
export SLEPC_DIR=/home/jess/develop/slepc/slepc
export PETSC_ARCH=arch-linux-c-debug
```
Переменная `PETSC_ARCH` задается в зависимости от вариантов сборки, которые доступны в пакете/репозитории. По умолчанию в репозитории есть `arch-linux-c-debug`.

3. Конфигурация, компиляция и проверка:
```
./configure 
make all
make check
```

### Установка SLEPc

4. Перейдите в корень репозитория SLEPc.

5. Конфигурация, компиляция и проверка:
```
./configure 
make all
make check
```

Если на каком-то этапе у вас что-то не получилось - попробуйте повторить шаг 2.

### Проверка установки

Запустите сборку тестовой программы, которая находится в одной папке с этим руководством:
```
make
```

## Windows

TODO (ставьте линукс)

### Полезные ссылки

* [SLEPc GitHub](https://github.com/slepc/slepc?tab=readme-ov-file)
* [PETSc download page](https://petsc.org/release/install/download/#recommended-obtain-release-version-with-git)
* [Installation guide](https://www.youtube.com/watch?v=hqVpaZ-axNk&t=428s)
* [Program compiling guide](https://www.youtube.com/watch?v=YeusG3o68vk)
