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

# Инструкция по установке PETSc и SLEPc и запуску программы на windows(Ubuntu (WSL))

## Введение
Эта инструкция поможет вам скачать, установить и настроить библиотеки **PETSc** и **SLEPc** в среде **Ubuntu/WSL**. Также будет показано, как скомпилировать и запустить программу на их основе.

## Шаг 1: Установка необходимых инструментов
Перед началом установки PETSc и SLEPc, убедитесь, что у вас установлены все необходимые утилиты для сборки:

```bash
sudo apt update
sudo apt install build-essential git python3 python3-pip libopenmpi-dev openmpi-bin
```

## Шаг 2: Скачивание PETSc
### 1. Скачайте PETSc
Скачайте PETSc с официального сайта или с помощью Git:

```bash
git clone https://gitlab.com/petsc/petsc.git
```

### 2. Перейдите в директорию с PETSc

```bash
cd petsc
```

### 3. Запустите конфигурацию PETSc
Вы можете скачать и установить все необходимые зависимости для PETSc (например, MPICH и BLAS/LAPACK) автоматически с помощью следующей команды:

```bash
./configure --download-fblaslapack --download-mpich
```

После конфигурации вам будет предложено собрать библиотеки.

### 4. Сборка PETSc

```bash
make all
```

### 5. Тестирование PETSc
Убедитесь, что PETSc работает корректно:

```bash
make test
```

## Шаг 3: Скачивание и установка SLEPc
### 1. Скачайте SLEPc
Для PETSc 3.18 используйте следующую команду:

```bash
git clone -b v3.18.3 https://gitlab.com/slepc/slepc.git
```

### 2. Перейдите в директорию SLEPc

```bash
cd slepc
```

### 3. Конфигурирование SLEPc
Убедитесь, что пути к PETSc правильно указаны через переменные окружения:

```bash
export PETSC_DIR=/path/to/petsc
export PETSC_ARCH=arch-linux-c-debug
```

Затем запустите конфигурацию:

```bash
./configure
```

### 4. Сборка SLEPc

```bash
make all
```

### 5. Тестирование SLEPc
Убедитесь, что SLEPc работает корректно:

```bash
make test
```

## Шаг 4: Написание программы с использованием PETSc и SLEPc
Создайте к примеру любой файл `main.c` и вставьте код.

## Шаг 5: Компиляция программы
Скомпилируйте программу с помощью `mpicc`:

```bash
mpicc -o eigenproblem main.c -I$PETSC_DIR/include -I$PETSC_DIR/$PETSC_ARCH/include -I$SLEPC_DIR/include -I$SLEPC_DIR/$PETSC_ARCH/include -L$PETSC_DIR/$PETSC_ARCH/lib -lpetsc -L$SLEPC_DIR/$PETSC_ARCH/lib -lslepc
```

## Шаг 6: Запуск программы
Запустите программу с использованием `mpirun`:

```bash
mpirun -np 2 ./eigenproblem
```

## Шаг 7: Исправление ошибок с библиотеками
Если при запуске программы вы видите ошибку типа `libpetsc.so.3.18: cannot open shared object file`, выполните следующие действия:

### 1. Добавьте библиотеки в переменную `LD_LIBRARY_PATH`

```bash
export LD_LIBRARY_PATH=$PETSC_DIR/$PETSC_ARCH/lib:$SLEPC_DIR/$PETSC_ARCH/lib:$LD_LIBRARY_PATH
```

### 2. Запустите программу снова

```bash
mpirun -np 2 ./eigenproblem
```

Теперь программа должна корректно запускаться.

## Заключение
Эта инструкция поможет вам настроить PETSc и SLEPc и выполнить вашу программу. Если возникнут дополнительные ошибки, проверьте правильность путей к библиотекам и настройки окружения.


### Полезные ссылки

* [SLEPc GitHub](https://github.com/slepc/slepc?tab=readme-ov-file)
* [PETSc download page](https://petsc.org/release/install/download/#recommended-obtain-release-version-with-git)
* [Installation guide](https://www.youtube.com/watch?v=hqVpaZ-axNk&t=428s)
* [Program compiling guide](https://www.youtube.com/watch?v=YeusG3o68vk)
