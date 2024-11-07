## Введение

Так как этот алгоритм в **slepc** — по сути, свёртка библиотеки **ScaLAPACK**, то для общих сведений кратко опишем саму библиотеку **ScaLAPACK**, исходя из руководства пользователя ([см. References[1]](https://netlib.org/scalapack/slug/node1.html#SECTION01000000000000000000)):

- **ScaLAPACK** — это библиотека высокопроизводительных процедур линейной алгебры для компьютеров MIMD с распределенной памятью и сетей рабочих станций, поддерживающих PVM и/или MPI.
- **Компоненты программного обеспечения ScaLAPACK**:
  - LAPACK
  - BLAS
  - PBLAS
  - BLACS
  - MPI

**NOTE**: Компоненты бывают разные: ([см. References[4]](https://netlib.org/scalapack/slug/node11.html#SECTION04130000000000000000)).

### LAPACK

- **ScaLAPACK** — это продолжение проекта LAPACK, в результате коего разработали и создали аналогичное (LAPACK) программное обеспечение для рабочих станций, векторных суперкомпьютеров и параллельных компьютеров с общей памятью — т.е. решили масштабировать LAPACK на большие архитектуры компьютеров.
- **Схожести ScaLAPACK и LAPACK**:
  - Обе библиотеки содержат процедуры для решения систем линейных уравнений, задач наименьших квадратов и задач на собственные значения.
  - Целями обоих проектов (ScaLAPACK и LAPACK) являются: ([см. References[2]](https://netlib.org/scalapack/slug/node9.html#SECTION04110000000000000000)).

**NOTE**: LAPACK будет работать на любой машине, где доступны BLAS, а ScaLAPACK будет работать на любой машине, где доступны и BLAS, и BLACS.

**NOTE2**: Процедуры в библиотеках хоть и схожи, но не идентичны.

### BLAS — обеспечивает работу LAPACK

- Это стандарт, используемый в данном контексте для LAPACK — в BLAS крайне оптимизированы основные операции линейной алгебры, а LAPACK, как больший пакет, пользуется качественными оптимизациями из BLAS.
- Важной целью BLAS является обеспечение уровня переносимости для вычислений.

### PBLAS — обеспечивает параллелизм

- Это реализация BLAS, обеспечивающая вычислительную основу уже для ScaLAPACK. ([References[6]](https://netlib.org/scalapack/slug/node14.html#SECTION04133000000000000000)).

### BLACS — обеспечивает специфическое общение между процессами

- Это библиотека передачи сообщений, разработанная для линейной алгебры.
- Простое умозаключение по ([References[7]](https://netlib.org/scalapack/slug/node15.html#SECTION04134000000000000000)): эта библиотека обеспечивает общение параллельных процессов — процессы «начинают общаться на языке линейной алгебры как на родном».

### MPI — обеспечивает общение между процессами ([References[8]](https://ru.wikipedia.org/wiki/Message_Passing_Interface))

---

## Описание алгоритма

По сути, алгоритм представляет собой реализацию свёртки ScaLAPACK.

См. описание алгоритмов от Коли. Причина — кажется, алгоритмы те же, смысл существования этого «алгоритма: ScaLAPACK» в slepc — в масштабируемости библиотеки.

---

## Имплементация

Покамест забиваем.

---

## References

1. [https://netlib.org/scalapack/slug/node1.html#SECTION01000000000000000000](https://netlib.org/scalapack/slug/node1.html#SECTION01000000000000000000) — user’s guide for Scalapack
2. [https://netlib.org/scalapack/slug/node9.html#SECTION04110000000000000000](https://netlib.org/scalapack/slug/node9.html#SECTION04110000000000000000) — user’s guide: Scalapack
3. [https://github.com/Reference-ScaLAPACK/scalapack](https://github.com/Reference-ScaLAPACK/scalapack) — репозиторий ScaLAPACK.
4. [https://netlib.org/scalapack/slug/node11.html#SECTION04130000000000000000](https://netlib.org/scalapack/slug/node11.html#SECTION04130000000000000000) — ScaLAPACK Software components.
5. [https://www.netlib.org/lapack/#_users_guide](https://www.netlib.org/lapack/#_users_guide) — страница lapack в репозитории netlib.
6. [https://netlib.org/scalapack/slug/node14.html#SECTION04133000000000000000](https://netlib.org/scalapack/slug/node14.html#SECTION04133000000000000000) — PBLAS
7. [https://netlib.org/scalapack/slug/node15.html#SECTION04134000000000000000](https://netlib.org/scalapack/slug/node15.html#SECTION04134000000000000000) — BLACS (Basic Linear Algebra Communication Subprograms)
8. [https://ru.wikipedia.org/wiki/Message_Passing_Interface](https://ru.wikipedia.org/wiki/Message_Passing_Interface) — MPI
