# Cross.c
# 1. Введение
Метод svd CrossMatrix основан на использовании кросс-матрицы:
C = A * $A^T$, где A - исходная матрица, $A^T$ - ее транспонированная матрица.


Из SLEPc использует собственные значения для определения структуры матрицы, вычисления сингулярных чисел и создания эффективной аппроксимации. Этот подход позволяет существенно сокращать размерность данных, сохраняя при этом основную информацию, что делает его полезным в различных областях науки и техники.

Задача решается преимущественно для гиперболического и обобщеннго типов. Что дает нам воспользоваться формулой:

$\sigma = \sqrt{|\lambda|}$, где $\lambda$ - собственное значение, $\sigma$ - сингялурное соотвествующее значение.

Матрица A является гиперболически типа, если выполняются следующие условия:
1. Все собственные значения λ_i действительные.
2. Среди этих собственных значений есть хотя бы одно положительное и хотя бы одно отрицательное.

Другими словами, если среди собственных значений матрицы есть как положительные, так и отрицательные, то такая матрица является гиперболическим типа.

Задача сингулярного разложения для пары матриц возникает в так называемой общей задаче сингулярных значений (Generalized Singular Value Decomposition, GSVD). В отличие от стандартного SVD, где рассматривается одна матрица, в GSVD рассматриваются две матрицы одновременно. Задача формулируется следующим образом:

Даны две матрицы A и B размера m×n и p×n, соответственно. Требуется найти такие унитарные матрицы U, V и диагональную матрицу Σ, что:

### $A = U * Σ_A * X^T$
### $B = V * Σ_B * X^T$
Где X и Y — унитарные матрицы, а Σ — диагональная матрица с сингулярными значениями.
#### Роль матриц A и B:
##### - Представляет первую и вторую матрицу, участвующую в задаче GSVD.
##### - Используются для формирования первого и второго набора сингулярных векторов.
#### Роль векторов u и v:
##### - Соответствует первому и второму набору сингулярных векторов, полученных из матрицы A и B.
##### - Являются результатом операции MatMult.

В библиотеке SLEPc реализованы для поиска собственных значений такие методы, как:
- Krylov-Schur Method (KSP). 

Этот метод является одним из основных в PETSc для решения систем линейных уравнений. Он использует итерационные методы, такие как GMRES, BiCGStab и другие, для нахождения собственных значений.

- Generalized Minimal Residual Method (GMRES) 

GMRES — это итерационный метод, используемый для решения систем линейных уравнений. Он эффективен для плохо обусловленных систем.

- Biconjugate Gradient Stabilized Method (BiCGStab)

BiCGStab — это еще один итерационный метод, который используется для решения системы линейных уравнений. Он также эффективен для плохо обусловленных систем.

- Conjugate Gradient Method (CG)

CG-метод используется для решения систем линейных уравнений с симметричными положительно определенными матрицами.

- Preconditioned Conjugate Gradient Method (PCG)

PCG-метод — это модификация CG-метода, использующая предобусловливание для улучшения сходимости.

##### Выбор конкретного метода зависит от специфики вашей задачи и требований к производительности и точности.

Параметр swapped в SLEPc отвечает за состояние, когда решение задачи сингулярных значений было получено таким образом, что матрицы были решены в порядке, отличающемся от того, какой порядок ожидался изначально. Это может произойти, например, если решение было найдено с использованием метода, который решает задачу для одного порядка матриц, тогда как ожидается другой порядок.

# 2. Описание алгоритма. 
## Здесь мы рассмотрим 2 алгоритма-функции: 

## Функция SVDSolve_Cross - Решает проблему сингулярного значения. 
##### - Вначале находится собственное значение для кросс матрицы одним из выше представленных методов. Так же получаем количество сходящихся собственных значений (nconv), количество выполненных итераций (its) и причину сходимости (reason).
##### - Затем проходясь по каждому собственному числу. Мы берем только вещественную часть собственного числа.
##### - Определяем принадлежит ли наша задача к гиперболическому типу, если да то $\sigma = \sqrt{|\sigma|}$. 
##### - Если нет, то проверяем проверяется больше ли $\sigma>-10* \epsilon$, где $\epsilon$ машинная точность. И выдается придупреждение об этом.
##### - Затем если $\sigma<0.0$ выполняется неравенство, то $\sigma = 0.0$
##### - Далее, берется корень из вещественно части $\sigma = \sqrt{\sigma}$

```
static PetscErrorCode SVDSolve_Cross(SVD svd)
{
  SVD_CROSS      *cross = (SVD_CROSS*)svd->data;
  PetscInt       i;
  PetscScalar    lambda;
  PetscReal      sigma;
  PetscFunctionBegin;
  PetscCall(EPSSolve(cross->eps));
  PetscCall(EPSGetConverged(cross->eps,&svd->nconv));
  PetscCall(EPSGetIterationNumber(cross->eps,&svd->its));
  PetscCall(EPSGetConvergedReason(cross->eps,(EPSConvergedReason*)&svd->reason));
  for (i=0;i<svd->nconv;i++) {
    PetscCall(EPSGetEigenvalue(cross->eps,i,&lambda,NULL));
    sigma = PetscRealPart(lambda);
    if (svd->ishyperbolic) svd->sigma[i] = PetscSqrtReal(PetscAbsReal(sigma));
    else {
      PetscCheck(sigma>-10*PETSC_MACHINE_EPSILON,PetscObjectComm((PetscObject)svd),PETSC_ERR_FP,"Negative eigenvalue computed by EPS: %g",(double)sigma);
      if (sigma<0.0) {
        PetscCall(PetscInfo(svd,"Negative eigenvalue computed by EPS: %g, resetting to 0\n",(double)sigma));
        sigma = 0.0;
      }
      svd->sigma[i] = PetscSqrtReal(sigma);
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}
```

# Функция SVDComputeVectors_Cross(SVD svd) - вычисление сингулярных векторов для обобщённых задач, гиперболические задачи и стандартные задачи, возникающие в рамках метода сингулярного разложения.

#### 1. Если задача ставится общая задача, то:
##### - создаются 2 вектора нужного размера под матирцы A, B
##### - Для каждого собственного значения: Извлекается i-й столбец из матрицы V, сохраняется в x и соответствующее ему собственное значение lambda. Затем выполняем равенства u = A * x. v = B * x.
##### - Нормализуем полученные векторы u, v.
##### - Вычисляем Скаляр alpha на основе отношения между собственными значениями и нормами векторов.
##### - Собственный вектор x масштабируется на величину alpha.
##### - Вектор x возвращается обратно в матрицу V
##### - Далее формируется вектор uv = (u,v), который состоит из полученных собственных векторов u, v на предыдущих шагах. Освобождается память от временных переменных.
##### 2. Если задача гиперболического типа, Swapped = true, 
##### - то получаем оператор Omega из EPS, создается вектор w подходящего размера для матрицы Omega, новый вектор для хранения сигнатуры.
##### - Для каждого собственного вектора: Извлечение i-го столбца из матрицы V и сохранение в него собственный вектор. 
##### - Умножение матрицы Omega на вектор v, результат сохраняется в w
##### - Вычисление скалярного произведения векторов v и w, результат сохраняем в alpha.
##### - Вычисляется сигнатура от вещественной части alpha.
##### - Затем $\alpha = \frac{1}{\sqrt{|\alpha|}}$
##### - Вектор w масштабируется с коэффициентом alpha, чтобы сохранить правильную норму. И копируется в вектор v.
##### - Значение сигнатуры сохраняется в векторе omega2, а индекс вектора в массиве varray обновляется.
##### - Продолжается вычисление оставшихся векторов с помощью функции SVDComputeVectors_Left
##### 3. Если задача негиперболического типа, Swapped = false,
##### - Для каждого собственного вектора: Извлекается i-й столбец из матрицы V и сохраняется в v.
##### - Извлекается i-й собственный вектор из объекта EPS и записывается в v.
##### - Возвращение вектора v обратно в матрицу V
##### - Вызов функции для вычисления левых сингулярных векторов.
```
static PetscErrorCode SVDComputeVectors_Cross(SVD svd)
{
  SVD_CROSS         *cross = (SVD_CROSS*)svd->data;
  PetscInt          i,mloc,ploc;
  Vec               u,v,x,uv,w,omega2=NULL;
  Mat               Omega;
  PetscScalar       *dst,alpha,lambda,*varray;
  const PetscScalar *src;
  PetscReal         nrm;
  PetscFunctionBegin;
  if (svd->isgeneralized) {
    PetscCall(MatCreateVecs(svd->A,NULL,&u));
    PetscCall(VecGetLocalSize(u,&mloc));
    PetscCall(MatCreateVecs(svd->B,NULL,&v));
    PetscCall(VecGetLocalSize(v,&ploc));
    for (i=0;i<svd->nconv;i++) {
      PetscCall(BVGetColumn(svd->V,i,&x));
      PetscCall(EPSGetEigenpair(cross->eps,i,&lambda,NULL,x,NULL));
      PetscCall(MatMult(svd->A,x,u));     /* u_i*c_i/alpha = A*x_i */
      PetscCall(VecNormalize(u,NULL));
      PetscCall(MatMult(svd->B,x,v));     /* v_i*s_i/alpha = B*x_i */
      PetscCall(VecNormalize(v,&nrm));    /* ||v||_2 = s_i/alpha   */
      alpha = 1.0/(PetscSqrtReal(1.0+PetscRealPart(lambda))*nrm);    /* alpha=s_i/||v||_2 */
      PetscCall(VecScale(x,alpha));
      PetscCall(BVRestoreColumn(svd->V,i,&x));
      /* copy [u;v] to U[i] */
      PetscCall(BVGetColumn(svd->U,i,&uv));
      PetscCall(VecGetArrayWrite(uv,&dst));
      PetscCall(VecGetArrayRead(u,&src));
      PetscCall(PetscArraycpy(dst,src,mloc));
      PetscCall(VecRestoreArrayRead(u,&src));
      PetscCall(VecGetArrayRead(v,&src));
      PetscCall(PetscArraycpy(dst+mloc,src,ploc));
      PetscCall(VecRestoreArrayRead(v,&src));
      PetscCall(VecRestoreArrayWrite(uv,&dst));
      PetscCall(BVRestoreColumn(svd->U,i,&uv));
    }
    PetscCall(VecDestroy(&v));
    PetscCall(VecDestroy(&u));
  } else if (svd->ishyperbolic && svd->swapped) {  /* was solved as GHIEP, set u=Omega*u and normalize */
    PetscCall(EPSGetOperators(cross->eps,NULL,&Omega));
    PetscCall(MatCreateVecs(Omega,&w,NULL));
    PetscCall(VecCreateSeq(PETSC_COMM_SELF,svd->ncv,&omega2));
    PetscCall(VecGetArrayWrite(omega2,&varray));
    for (i=0;i<svd->nconv;i++) {
      PetscCall(BVGetColumn(svd->V,i,&v));
      PetscCall(EPSGetEigenvector(cross->eps,i,v,NULL));
      PetscCall(MatMult(Omega,v,w));
      PetscCall(VecDot(v,w,&alpha));
      svd->sign[i] = PetscSign(PetscRealPart(alpha));
      varray[i] = svd->sign[i];
      alpha = 1.0/PetscSqrtScalar(PetscAbsScalar(alpha));
      PetscCall(VecScale(w,alpha));
      PetscCall(VecCopy(w,v));
      PetscCall(BVRestoreColumn(svd->V,i,&v));
    }
    PetscCall(BVSetSignature(svd->V,omega2));
    PetscCall(VecRestoreArrayWrite(omega2,&varray));
    PetscCall(VecDestroy(&omega2));
    PetscCall(VecDestroy(&w));
    PetscCall(SVDComputeVectors_Left(svd));
  } else {
    for (i=0;i<svd->nconv;i++) {
      PetscCall(BVGetColumn(svd->V,i,&v));
      PetscCall(EPSGetEigenvector(cross->eps,i,v,NULL));
      PetscCall(BVRestoreColumn(svd->V,i,&v));
    }
    PetscCall(SVDComputeVectors_Left(svd));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}
```

# 3. Имплементация.
### Импортится обработчик ошибок:
``` #include <slepc/private/svdimpl.h> ``` 

### Структура SVD_CROSS [1]
```
typedef struct {
  PetscBool explicitmatrix;
  EPS       eps;
  PetscBool usereps;
  Mat       C,D;
} SVD_CROSS;
```
Которая состоит: 
 - explicitmatrix - параметр, который отвечает является ли матрица явной или неявной, т.е. задана она в явном виде или нет. Что такое неявная матрица дальше.
 - eps - Решатель задач на собственные значения (EPS) — это основной объект, предоставляемый slepc. Он используется для задания линейной задачи на собственные значения в стандартной или обобщенной форме и обеспечивает единый и эффективный доступ ко всем линейным решателям собственных значений, включенным в пакет. Концептуально уровень абстракции, занимаемый EPS, аналогичен другим решателям в petsc, таким как KSP, для решения линейных систем уравнений. Модуль EPS в slepc используется аналогично модулям petsc, таким как KSP. Вся информация, связанная с задачей на собственные значения, обрабатывается через контекстную переменную. Доступны обычные функции управления объектами (EPSCreate, EPSDestroy, EPSView, EPSSetFromOptions). Кроме того, объект EPS предоставляет функции для установки нескольких параметров, таких как количество вычисляемых собственных значений, размерность подпространства, часть интересующего спектра, запрашиваемый допуск или максимальное количество разрешенных итераций. [1]

- usereps - переменная типа bool, для настройки индивидуальных параметров какого-то алгоритма.
- C,D - какие-то матрицы.

Примеры неявных матриц: 
1. Разреженные матрицы
Многие задачи требуют работы с разреженными матрицами, у которых большинство элементов равны нулю. Вместо того чтобы хранить всю матрицу целиком, используют специальные структуры данных, такие как списки смежностей или формат COO (координатный формат), которые хранят только ненулевые элементы вместе с их индексами. Таким образом, матрица становится неявной, поскольку большая часть её данных не хранится явно.

2. Квантование матриц
Некоторые методы квантования представляют матрицы в сжатом виде, используя ограниченное количество бит для каждого элемента. Это позволяет уменьшить размер матрицы за счёт точности, и в результате мы получаем неявное представление матрицы.

3. Алгоритмическое задание матриц
Иногда матрицы могут задаваться алгоритмами, которые генерируют элементы по мере необходимости. Например, если матрица является результатом некоторого процесса, её элементы могут вычисляться динамически, а не храниться в памяти заранее.

4. Использование функций для генерации элементов
В некоторых случаях элементы матрицы можно задать функцией, которая возвращает значение элемента по его индексу. Такой подход особенно полезен, когда необходимо работать с очень большими матрицами или с бесконечными последовательностями.

Преимущества неявного представления матриц:
 - Экономия памяти: Для больших матриц хранение всех элементов может потребовать значительных ресурсов, тогда как неявная форма позволяет избежать этого.
 - Более эффективное выполнение операций: Некоторые операции над матрицами могут выполняться быстрее, если они реализованы через функции, а не через прямое обращение к каждому элементу.
 - Упрощение логики программы: В некоторых случаях работа с неявными матрицами упрощает код и делает его более понятным.

Недостатки:
 - Ограничения на некоторые виды операций: Не все операции можно эффективно выполнять с неявно заданными матрицами. Например, сложение двух неявно определённых матриц может требовать их полного раскрытия.
 - Сложность реализации: Работа с неявным представлением требует дополнительных усилий при разработке и тестировании алгоритмов.

### Структура SVD_CROSS_SHELL[1]
```
typedef struct {
  Mat       A,AT;
  Vec       w,diag,omega;
  PetscBool swapped;
} SVD_CROSS_SHELL;
```
Которая состоит из: 
- A - проблемная матрица
- AT - транспонированная матрица;
- w -- вектор, который используется для сохранения промежуточных результатов вычислений в процессе декомпозиции задачи сингулярных значений.,
- diag - диагональные элементы матрицы,
- omega - сигнатура для гиперболических задач;
- swapped - U и V поменялись местами (M<N). U,V - левый и правый сингулярные векторы;

# PetscErrorCode MatMult_Cross - функция которая принимает на вход примает входную матрицу, входной вектор, выходной вектор. Которая позволяет умножить матрицу B на вектор x.

```
static  PetscErrorCode MatMult_Cross(Mat B,Vec x,Vec y)
{
  SVD_CROSS_SHELL *ctx; 
  
  PetscFunctionBegin;
  PetscCall(MatShellGetContext(B,&ctx));
  PetscCall(MatMult(ctx->A,x,ctx->w));
  if (ctx->omega && !ctx->swapped) PetscCall(VecPointwiseMult(ctx->w,ctx->w,ctx->omega));
  PetscCall(MatMult(ctx->AT,ctx->w,y));
  PetscFunctionReturn(PETSC_SUCCESS);
}
```
 
# Функция MatGetDiagonal_Cross, которая вычисляет диагональ матрицы и сохраняет её в вектор. Mat B - входная матрица, выходной вектор d (диагональ матрицы)
```
static PetscErrorCode MatGetDiagonal_Cross(Mat B,Vec d)
{
  SVD_CROSS_SHELL   *ctx;
  PetscMPIInt       len;
  PetscInt          N,n,i,j,start,end,ncols;
  PetscScalar       *work1,*work2,*diag;
  const PetscInt    *cols;
  const PetscScalar *vals;

  PetscFunctionBegin;
  PetscCall(MatShellGetContext(B,&ctx));
  if (!ctx->diag) {
    /* compute diagonal from rows and store in ctx->diag */
    PetscCall(VecDuplicate(d,&ctx->diag));
    PetscCall(MatGetSize(ctx->A,NULL,&N));
    PetscCall(MatGetLocalSize(ctx->A,NULL,&n));
    PetscCall(PetscCalloc2(N,&work1,N,&work2));
    if (ctx->swapped) {
      PetscCall(MatGetOwnershipRange(ctx->AT,&start,&end));
      for (i=start;i<end;i++) {
        PetscCall(MatGetRow(ctx->AT,i,&ncols,NULL,&vals));
        for (j=0;j<ncols;j++) work1[i] += vals[j]*vals[j];
        PetscCall(MatRestoreRow(ctx->AT,i,&ncols,NULL,&vals));
      }
    } else {
      PetscCall(MatGetOwnershipRange(ctx->A,&start,&end));
      for (i=start;i<end;i++) {
        PetscCall(MatGetRow(ctx->A,i,&ncols,&cols,&vals));
        for (j=0;j<ncols;j++) work1[cols[j]] += vals[j]*vals[j];
        PetscCall(MatRestoreRow(ctx->A,i,&ncols,&cols,&vals));
      }
    }
    PetscCall(PetscMPIIntCast(N,&len));
    PetscCallMPI(MPIU_Allreduce(work1,work2,len,MPIU_SCALAR,MPIU_SUM,PetscObjectComm((PetscObject)B)));
    PetscCall(VecGetOwnershipRange(ctx->diag,&start,&end));
    PetscCall(VecGetArrayWrite(ctx->diag,&diag));
    for (i=start;i<end;i++) diag[i-start] = work2[i];
    PetscCall(VecRestoreArrayWrite(ctx->diag,&diag));
    PetscCall(PetscFree2(work1,work2));
  }
  PetscCall(VecCopy(ctx->diag,d));
  PetscFunctionReturn(PETSC_SUCCESS);
}
```
1. Объявление необходимых параметров: 
2.
```
PetscCall(MatShellGetContext(B,&ctx))
``` 
Возвращает предоставленный пользователем контекст, связанный с матрицей оболочки MATSHELL. Входной параметр: B - матрица, должна быть создана с помощью MatCreateShell() Выходной параметр: ctx - предоставленный пользователем контекст.
3.
```
if (!ctx->diag) 
``` 
Проверяем есть ли вычисленная диагональ матрицы в контексте. То выполняются пункт 4-6.
4.1
```
PetscCall(VecDuplicate(d,&ctx->diag));
``` 
 Создает новый вектор того же типа, что и существующий вектор.
4.2
```
    PetscCall(MatGetSize(ctx->A,NULL,&N));
``` 
Возвращает количество столбцов в матрице.
4.3 
```
PetscCall(MatGetLocalSize(ctx->A,NULL,&n)); 
``` 
Получает размерность матрицы ctx->A и сохраняет результат в переменной n
4.3 
```
PetscCall(PetscCalloc2(N,&work1,N,&work2)); 
```
Выделяет памяти для 2 блоков размера .
5) 
```
if (ctx->swapped)
``` 
выбор вычисления взавизимости от swapped.
5.1) 
```
PetscCall(MatGetOwnershipRange(ctx->AT,&start,&end));
``` 
Здесь берется транспонированная матрица. - Для матриц, владеющих значениями по строкам, исключая MATELEMENTAL и MATSCALAPACK, возвращает диапазон строк матрицы, принадлежащих этому процессу MPI.
```
      for (i=start;i<end;i++) {
```
```
PetscCall(MatGetRow(ctx->AT,i,&ncols,NULL,&vals));
``` 
Получает строку матрицы.
```
for (j=0;j<ncols;j++) work1[i] += vals[j]*vals[j];
``` 
Для каждой строки матрицы извлекаются ненулевые элементы и их индексы (в случае обычной матрицы), затем каждый элемент возводится в квадрат и добавляется к соответствующему индексу в рабочем массиве work1.
```
PetscCall(MatRestoreRow(ctx->AT,i,&ncols,NULL,&vals));
``` 
Освобождает любое временное пространство, выделенное MatGetRow()

5.2) Делается все тоже самое как в пункте 5.1) Только для обычной матрица A, а не транспонированной.
```
      PetscCall(MatGetOwnershipRange(ctx->A,&start,&end));
      for (i=start;i<end;i++) {
        PetscCall(MatGetRow(ctx->A,i,&ncols,&cols,&vals));
        for (j=0;j<ncols;j++) work1[cols[j]] += vals[j]*vals[j];
        PetscCall(MatRestoreRow(ctx->A,i,&ncols,&cols,&vals));
      }
```
6.1) 
``` 
PetscCall(PetscMPIIntCast(N,&len));
``` 
преобразует MPI_Count, PetscInt, PetscCount или PetscInt64 в PetscMPIInt (который всегда имеет размер 32 бита), генерирует ошибку, если PetscMPIInt недостаточно велик для хранения числа.
6.2) 
```
PetscCallMPI(MPIU_Allreduce(work1,work2,len,MPIU_SCALAR,MPIU_SUM,PetscObjectComm((PetscObject)B)));
``` 
Замена для MPI_Allreduce(), которая (1) выполняет операции MPIU_INT с одним счетчиком в PetscInt64 для обнаружения переполнений целых чисел и (2) пытается определить, происходит ли вызов из всех рангов MPI в одном и том же месте в коде PETSc. Это помогает обнаружить ошибки, когда разные ранги MPI следуют разным путям кода, что приводит к непоследовательным и некорректным вызовам.
```
PetscMPIInt MPIU_Allreduce(void *indata,void *outdata,PetscCount count,MPI_Datatype dtype, MPI_Op op, MPI_Comm comm);
```
Входные параметры:
a - указатель на входные данные, которые необходимо сократить;
count - количество элементов данных MPI в a и b;
dtype - тип данных MPI, например MPI_INT;
op - операция MPI, например MPI_SUM;
comm - коммуникатор MPI, на котором выполняется операция.
Выходной параметр:
b - сокращенные значения;
6.3) 
```PetscCall(VecGetOwnershipRange(ctx->diag,&start,&end));
``` 
Возвращает диапазон индексов, принадлежащих этому процессу. Вектор размещается с первыми n1 элементами на первом процессоре, следующими n2 элементами на втором и т. д. Для некоторых параллельных макетов этот диапазон может быть неопределенным.
6.4)
```
PetscCall(VecGetArrayWrite(ctx->diag,&diag));
```
Возвращает указатель на непрерывный массив, который БУДЕТ содержать часть векторных данных этого процесса MPI.
6.5)
```
for (i=start;i<end;i++) diag[i-start] = work2[i];
PetscCall(VecRestoreArrayWrite(ctx->diag,&diag));
``` 
Диагональные элементы копируются из массива work2 в вектор ctx->diag.
6.7)
```
PetscCall(PetscFree2(work1,work2));
``` 
Освобождается выделенная ранее память для рабочих массивов.
7) 
```
PetscCall(VecCopy(ctx->diag,d));
PetscFunctionReturn(PETSC_SUCCESS);
```
Делается копия вектора.

# MatDestroy_Cross(Mat B) - Деструктор матрицы B.
```
static PetscErrorCode MatDestroy_Cross(Mat B)
{
  SVD_CROSS_SHELL *ctx;

  PetscFunctionBegin;
  PetscCall(MatShellGetContext(B,&ctx));
  PetscCall(VecDestroy(&ctx->w));
  PetscCall(VecDestroy(&ctx->diag));
  PetscCall(PetscFree(ctx));
  PetscFunctionReturn(PETSC_SUCCESS);
}
```

# SVDSetUp_Cross(SVD svd) - Устанавливает все внутренние структуры данных, необходимые для выполнения решателя сингулярных значений.
``` 
static PetscErrorCode SVDSetUp_Cross(SVD svd)
{
  SVD_CROSS      *cross = (SVD_CROSS*)svd->data;
  ST             st;
  PetscBool      trackall,issinv,isks;
  EPSProblemType ptype;
  EPSWhich       which;
  Mat            Omega;
  MatType        Atype;
  PetscInt       n,N;

  PetscFunctionBegin;
  if (!cross->eps) PetscCall(SVDCrossGetEPS(svd,&cross->eps));
  PetscCall(MatDestroy(&cross->C));
  PetscCall(MatDestroy(&cross->D));
  PetscCall(SVDCrossGetProductMat(svd,svd->A,svd->AT,&cross->C));
  if (svd->isgeneralized) {
    PetscCall(SVDCrossGetProductMat(svd,svd->B,svd->BT,&cross->D));
    PetscCall(EPSSetOperators(cross->eps,cross->C,cross->D));
    PetscCall(EPSGetProblemType(cross->eps,&ptype));
    if (!ptype) PetscCall(EPSSetProblemType(cross->eps,EPS_GHEP));
  } else if (svd->ishyperbolic && svd->swapped) {
    PetscCall(MatGetType(svd->OP,&Atype));
    PetscCall(MatGetSize(svd->A,NULL,&N));
    PetscCall(MatGetLocalSize(svd->A,NULL,&n));
    PetscCall(MatCreate(PetscObjectComm((PetscObject)svd),&Omega));
    PetscCall(MatSetSizes(Omega,n,n,N,N));
    PetscCall(MatSetType(Omega,Atype));
    PetscCall(MatDiagonalSet(Omega,svd->omega,INSERT_VALUES));
    PetscCall(EPSSetOperators(cross->eps,cross->C,Omega));
    PetscCall(EPSSetProblemType(cross->eps,EPS_GHIEP));
    PetscCall(MatDestroy(&Omega));
  } else {
    PetscCall(EPSSetOperators(cross->eps,cross->C,NULL));
    PetscCall(EPSSetProblemType(cross->eps,EPS_HEP));
  }
  if (!cross->usereps) {
    PetscCall(EPSGetST(cross->eps,&st));
    PetscCall(PetscObjectTypeCompare((PetscObject)st,STSINVERT,&issinv));
    PetscCall(PetscObjectTypeCompare((PetscObject)cross->eps,EPSKRYLOVSCHUR,&isks));
    if (svd->isgeneralized && svd->which==SVD_SMALLEST) {
      if (cross->explicitmatrix && isks && !issinv) {  /* default to shift-and-invert */
        PetscCall(STSetType(st,STSINVERT));
        PetscCall(EPSSetTarget(cross->eps,0.0));
        which = EPS_TARGET_REAL;
      } else which = issinv?EPS_TARGET_REAL:EPS_SMALLEST_REAL;
    } else {
      if (issinv) which = EPS_TARGET_MAGNITUDE;
      else if (svd->ishyperbolic) which = svd->which==SVD_LARGEST?EPS_LARGEST_MAGNITUDE:EPS_SMALLEST_MAGNITUDE;
      else which = svd->which==SVD_LARGEST?EPS_LARGEST_REAL:EPS_SMALLEST_REAL;
    }
    PetscCall(EPSSetWhichEigenpairs(cross->eps,which));
    PetscCall(EPSSetDimensions(cross->eps,svd->nsv,svd->ncv,svd->mpd));
    PetscCall(EPSSetTolerances(cross->eps,svd->tol==(PetscReal)PETSC_DETERMINE?SLEPC_DEFAULT_TOL/10.0:svd->tol,svd->max_it));
    switch (svd->conv) {
    case SVD_CONV_ABS:
      PetscCall(EPSSetConvergenceTest(cross->eps,EPS_CONV_ABS));break;
    case SVD_CONV_REL:
      PetscCall(EPSSetConvergenceTest(cross->eps,EPS_CONV_REL));break;
    case SVD_CONV_NORM:
      if (svd->isgeneralized) {
        if (!svd->nrma) PetscCall(MatNorm(svd->OP,NORM_INFINITY,&svd->nrma));
        if (!svd->nrmb) PetscCall(MatNorm(svd->OPb,NORM_INFINITY,&svd->nrmb));
        PetscCall(EPSSetConvergenceTestFunction(cross->eps,EPSConv_Cross,svd,NULL));
      } else {
        PetscCall(EPSSetConvergenceTest(cross->eps,EPS_CONV_NORM));break;
      }
      break;
    case SVD_CONV_MAXIT:
      SETERRQ(PetscObjectComm((PetscObject)svd),PETSC_ERR_SUP,"Maxit convergence test not supported in this solver");
    case SVD_CONV_USER:
      SETERRQ(PetscObjectComm((PetscObject)svd),PETSC_ERR_SUP,"User-defined convergence test not supported in this solver");
    }
  }
  SVDCheckUnsupported(svd,SVD_FEATURE_STOPPING);
  /* Transfer the trackall option from svd to eps */
  PetscCall(SVDGetTrackAll(svd,&trackall));
  PetscCall(EPSSetTrackAll(cross->eps,trackall));
  /* Transfer the initial space from svd to eps */
  if (svd->nini<0) {
    PetscCall(EPSSetInitialSpace(cross->eps,-svd->nini,svd->IS));
    PetscCall(SlepcBasisDestroy_Private(&svd->nini,&svd->IS));
  }
  PetscCall(EPSSetUp(cross->eps));
  PetscCall(EPSGetDimensions(cross->eps,NULL,&svd->ncv,&svd->mpd));
  PetscCall(EPSGetTolerances(cross->eps,NULL,&svd->max_it));
  if (svd->tol==(PetscReal)PETSC_DETERMINE) svd->tol = SLEPC_DEFAULT_TOL;

  svd->leftbasis = PETSC_FALSE;
  PetscCall(SVDAllocateSolution(svd,0));
  PetscFunctionReturn(PETSC_SUCCESS);
}
```
1) Обьявление нужных переменных
2) 
```
if (!cross->eps) PetscCall(SVDCrossGetEPS(svd,&cross->eps));
```
Проверяется наличие объекта EPS в структуре cross. Если его нет, он создается с помощью SVDCrossGetEPS.
```
PetscCall(MatDestroy(&cross->C));
PetscCall(MatDestroy(&cross->D));
```
Уничтожение подматриц C и D относящиеся к cross.
```
PetscCall(SVDCrossGetProductMat(svd,svd->A,svd->AT,&cross->C));
```
Объявление матрице cross->С, что она равна произведение матрицы A на ее транспонированную матрицу AT.
3) 
```
if (svd->isgeneralized) {
    PetscCall(SVDCrossGetProductMat(svd,svd->B,svd->BT,&cross->D));
    PetscCall(EPSSetOperators(cross->eps,cross->C,cross->D));
    PetscCall(EPSGetProblemType(cross->eps,&ptype));
    if (!ptype) PetscCall(EPSSetProblemType(cross->eps,EPS_GHEP));
```
Если задача является обобщенной, то 
```
PetscCall(SVDCrossGetProductMat(svd,svd->B,svd->BT,&cross->D));
```
Матрица cross->D равна произведению матрицы B на ее транспонированную матрицу BT.

```
 PetscCall(EPSSetOperators(cross->eps,cross->C,cross->D));
 ```
Эта функция устанавливает операторы для задачи собственных значений, которую решает объект типа EPS.
Параметры:
cross->eps: Объект типа EPS, для которого устанавливаются операторы.
cross->C: Матрица оператора, используемая для левой части уравнения собственных значений.
cross->D: Матрица оператора, используемая для правой части уравнения собственных значений (для обобщённых задач собственных значений)

```
PetscCall(EPSGetProblemType(cross->eps,&ptype));
```
Определение типа проблемы.
```
if (!ptype) PetscCall(EPSSetProblemType(cross->eps,EPS_GHEP));
```
Если тип проблемы не был задан, устанавливается тип EPS_GHEP (обобщенное собственное значение проблема).
4) 
```
else if (svd->ishyperbolic && svd->swapped) {
    PetscCall(MatGetType(svd->OP,&Atype));
    PetscCall(MatGetSize(svd->A,NULL,&N));
    PetscCall(MatGetLocalSize(svd->A,NULL,&n));
    PetscCall(MatCreate(PetscObjectComm((PetscObject)svd),&Omega));
    PetscCall(MatSetSizes(Omega,n,n,N,N));
    PetscCall(MatSetType(Omega,Atype));
    PetscCall(MatDiagonalSet(Omega,svd->omega,INSERT_VALUES));
    PetscCall(EPSSetOperators(cross->eps,cross->C,Omega));
    PetscCall(EPSSetProblemType(cross->eps,EPS_GHIEP));
    PetscCall(MatDestroy(&Omega));
  } else {
    PetscCall(EPSSetOperators(cross->eps,cross->C,NULL));
    PetscCall(EPSSetProblemType(cross->eps,EPS_HEP));
  }
```
Если задача является гиперболической и флаг swapped установлен. То вначале определяются параметры матрицы и тип оператора. Создается дополнительная матрица Omega, (аналогичная по размеру и типу оператора матрице А) заполняемая значениями вектора omega по диагонали. Эта матрица используется вместе с C для установки операторов в задаче собственных значений. Устанавливается тип задачи EPS_GHIEP (гиперболическая собственная величина проблема). По завершении работы матрица Omega уничтожается.
Если ни одно из предыдущих условий не выполнено, устанавливается стандартная задача собственных значений (EPS_HEP) с оператором C.
5) 
```
if (!cross->usereps) {
```
Условие проверяет, установлено ли значение usereps в структуре cross
```
    PetscCall(EPSGetST(cross->eps,&st));
```
EPSGetST определяет тип спектрального преобразования (ST) из объекта cross->eps и сохраняет её в переменной st.
```
    PetscCall(PetscObjectTypeCompare((PetscObject)st,STSINVERT,&issinv));
```
PetscObjectTypeCompare сравнивает тип объекта st с типом STSINVERT и присваивает результат сравнения переменной issinv.
```
    PetscCall(PetscObjectTypeCompare((PetscObject)cross->eps,EPSKRYLOVSCHUR,&isks));
```
Проверяет, является ли объект cross->eps типа EPSKRYLOVSCHUR, и результат сохранения в переменную isks
```
    if (svd->isgeneralized && svd->which==SVD_SMALLEST) {
```
Если задача является обобщённой (isgeneralized == true) и требуется найти наименьшие собственные значения (which == SVD_SMALLEST),
```
      if (cross->explicitmatrix && isks && !issinv) {  
```
Проверяется, явная ли матрица явной и используется ли метод Крылова-Шура (isks), и при этом не используется инвертированное преобразование (!issinv)
```
        PetscCall(STSetType(st,STSINVERT));
        PetscCall(EPSSetTarget(cross->eps,0.0));
        which = EPS_TARGET_REAL;
```
Если да, и при этом не используется инвертированное преобразование (!issinv), то устанавливается тип инвертированного преобразования (STSINVERT) и целевое значение для поиска собственных значений (EPSSetTarget).
```
      } else which = issinv?EPS_TARGET_REAL:EPS_SMALLEST_REAL;
    } else {
      if (issinv) which = EPS_TARGET_MAGNITUDE;
      else if (svd->ishyperbolic) which = svd->which==SVD_LARGEST?EPS_LARGEST_MAGNITUDE:EPS_SMALLEST_MAGNITUDE;
      else which = svd->which==SVD_LARGEST?EPS_LARGEST_REAL:EPS_SMALLEST_REAL;
    }
```
Если условие выше не выполняется, то проводится дальнейшая проверка:
Если используется инвертированное преобразование (issinv), то ищутся максимальные по модулю собственные значения (EPS_TARGET_MAGNITUDE).
Иначе, если задача гиперболическая (ishyperbolic), то ищут наибольшие или наименьшие по модулю собственные значения в зависимости от значения which.
В остальных случаях ищут наибольшие (EPS_LARGEST_REAL) или наименьшие (EPS_SMALLEST_REAL) действительные собственные значения.
```
    PetscCall(EPSSetWhichEigenpairs(cross->eps,which));
```
Устанавливает, какие именно собственные значения будут найдены (наибольшие, наименьшие и т.п.).
```
    PetscCall(EPSSetDimensions(cross->eps,svd->nsv,svd->ncv,svd->mpd));
```
Устанавливает количество требуемых собственных значений (nsv), количество начальных векторов (ncv) и максимальное число дополнительных векторов (mpd).
```
    PetscCall(EPSSetTolerances(cross->eps,svd->tol==(PetscReal)PETSC_DETERMINE?SLEPC_DEFAULT_TOL/10.0:svd->tol,svd->max_it));
```
Устанавливает допустимую погрешность и максимальное количество итераций для нахождения собственных значений.
```
    switch (svd->conv) {
    case SVD_CONV_ABS:
      PetscCall(EPSSetConvergenceTest(cross->eps,EPS_CONV_ABS));break;
    case SVD_CONV_REL:
      PetscCall(EPSSetConvergenceTest(cross->eps,EPS_CONV_REL));break;
    case SVD_CONV_NORM:
      if (svd->isgeneralized) {
        if (!svd->nrma) PetscCall(MatNorm(svd->OP,NORM_INFINITY,&svd->nrma));
        if (!svd->nrmb) PetscCall(MatNorm(svd->OPb,NORM_INFINITY,&svd->nrmb));
        PetscCall(EPSSetConvergenceTestFunction(cross->eps,EPSConv_Cross,svd,NULL));
      } else {
        PetscCall(EPSSetConvergenceTest(cross->eps,EPS_CONV_NORM));break;
      }
      break;
    case SVD_CONV_MAXIT:
      SETERRQ(PetscObjectComm((PetscObject)svd),PETSC_ERR_SUP,"Maxit convergence test not supported in this solver");
    case SVD_CONV_USER:
      SETERRQ(PetscObjectComm((PetscObject)svd),PETSC_ERR_SUP,"User-defined convergence test not supported in this solver");
    }
  }
```
В зависимости от установленного критерия сходимости (conv), выбираются разные методы проверки сходимости:
SVD_CONV_ABS: Абсолютная сходимость.
SVD_CONV_REL: Относительная сходимость.
SVD_CONV_NORM: Нормированная сходимость. Если задача обобщённая, дополнительно вычисляются нормы матриц OP и OPb.
Для SVD_CONV_MAXIT и SVD_CONV_USER выдаются сообщения об ошибке, поскольку данные критерии сходимости не поддерживаются этим решателем.
6) 
```
SVDCheckUnsupported(svd,SVD_FEATURE_STOPPING);
```
Перед настройкой объекта EPS проверяются неподдерживаемые параметры.
```
  PetscCall(SVDGetTrackAll(svd,&trackall));
```
Перенести опцию trackall из svd в eps
```
  PetscCall(EPSSetTrackAll(cross->eps,trackall));
```
Перенести перенести пространство из svd в eps

Установка начального пространства: 
```
  if (svd->nini<0) {
    PetscCall(EPSSetInitialSpace(cross->eps,-svd->nini,svd->IS));
    PetscCall(SlepcBasisDestroy_Private(&svd->nini,&svd->IS));
  }
```
Если svd->nini меньше нуля, это означает, что у нас есть заранее определенное начальное пространство для задачи собственных значений.
Функция EPSSetInitialSpace передает это начальное пространство объекту EPS для использования в процессе решения.
После передачи начального пространства, оно уничтожается с помощью вызова SlepcBasisDestroy_Private, чтобы избежать утечек памяти.
```
  PetscCall(EPSSetUp(cross->eps));
```
Вызов EPSSetUp завершает настройку объекта EPS, используя все ранее установленные параметры. После этого объект готов к решению задачи собственных значений.
```
  PetscCall(EPSGetDimensions(cross->eps,NULL,&svd->ncv,&svd->mpd));
```
Обновляет значения ncv (количество начальных векторов) и mpd (максимальное количество дополнительных векторов) в соответствии с тем, что было определено в ходе настройки объекта EPS.
```
  PetscCall(EPSGetTolerances(cross->eps,NULL,&svd->max_it));
```
Обновляет значение max_it (максимальное количество итераций) для объекта EPS.
```
  if (svd->tol==(PetscReal)PETSC_DETERMINE) svd->tol = SLEPC_DEFAULT_TOL;
```
Если значение svd->tol было установлено как PETSC_DETERMINE, то оно заменяется на значение по умолчанию SLEPC_DEFAULT_TOL
```
  svd->leftbasis = PETSC_FALSE;
```
Поле leftbasis объекта svd устанавливается в PETSC_FALSE, что говорит о том, что левое базисное решение не нужно сохранять.
```
  PetscCall(SVDAllocateSolution(svd,0));
```
Вызов SVDAllocateSolution выделяет память для хранения решения задачи собственных значений.
  PetscFunctionReturn(PETSC_SUCCESS);

# Фукнция EPSMonitor_Cross(EPS eps,PetscInt its,PetscInt nconv,PetscScalar *eigr,PetscScalar *eigi,PetscReal *errest,PetscInt nest,void *ctx) - контролирует процесс решения задачи собственных значений путем мониторинга ошибок на каждом уровне вложенности и обновления соответствующих метрик, таких как точность и ошибка.
```
static PetscErrorCode EPSMonitor_Cross(EPS eps,PetscInt its,PetscInt nconv,PetscScalar *eigr,PetscScalar *eigi,PetscReal *errest,PetscInt nest,void *ctx)
{
  PetscInt       i;
  SVD            svd = (SVD)ctx;
  PetscScalar    er,ei;
  ST             st;

  PetscFunctionBegin;
  PetscCall(EPSGetST(eps,&st));
  for (i=0;i<PetscMin(nest,svd->ncv);i++) {
    er = eigr[i]; ei = eigi[i];
    PetscCall(STBackTransform(st,1,&er,&ei));
    svd->sigma[i] = PetscSqrtReal(PetscAbsReal(PetscRealPart(er)));
    svd->errest[i] = errest[i];
  }
  PetscCall(SVDMonitor(svd,its,nconv,svd->sigma,svd->errest,nest));
  PetscFunctionReturn(PETSC_SUCCESS);
}
```
1) 
```
PetscCall(EPSGetST(eps,&st));
```
Вызывается функция EPSGetST, которая возвращает объект ST, отвечающий за тип преобразования, используемую в решении задачи собственных значений.
```
for (i=0;i<PetscMin(nest,svd->ncv);i++) {
    er = eigr[i]; ei = eigi[i];
```
Для каждого уровня вложенности выполняется обратное преобразование ошибок:

```
PetscCall(STBackTransform(st,1,&er,&ei));
```
Функция STBackTransform выполняет обратное преобразование ошибок от уровня вложения до исходного уровня.

```
svd->sigma[i] = PetscSqrtReal(PetscAbsReal(PetscRealPart(er)));
svd->errest[i] = errest[i];

```
Вектора sigma и errest обновляются согласно результатам преобразования ошибок
```
PetscCall(SVDMonitor(svd,its,nconv,svd->sigma,svd->errest,nest));
```
Вызывается функция SVDMonitor, которая осуществляет мониторинг решения задачи с учетом обновленных значений ошибок.

# SVDSetFromOptions_Cross(SVD svd,PetscOptionItems *PetscOptionsObject) - функция отвечает за настройку объекта SVD через интерфейс опций, позволяющий пользователям выбирать различные параметры для решения задачи сингулярного разложения методом cross SVD.
```
static PetscErrorCode SVDSetFromOptions_Cross(SVD svd,PetscOptionItems *PetscOptionsObject)
{
  PetscBool      set,val;
  SVD_CROSS      *cross = (SVD_CROSS*)svd->data;
  ST             st;

  PetscFunctionBegin;
  PetscOptionsHeadBegin(PetscOptionsObject,"SVD Cross Options");

    PetscCall(PetscOptionsBool("-svd_cross_explicitmatrix","Use cross explicit matrix","SVDCrossSetExplicitMatrix",cross->explicitmatrix,&val,&set));
    if (set) PetscCall(SVDCrossSetExplicitMatrix(svd,val));

  PetscOptionsHeadEnd();

  if (!cross->eps) PetscCall(SVDCrossGetEPS(svd,&cross->eps));
  if (!cross->explicitmatrix && !cross->usereps) {
    /* use as default an ST with shell matrix and Jacobi */
    PetscCall(EPSGetST(cross->eps,&st));
    PetscCall(STSetMatMode(st,ST_MATMODE_SHELL));
  }
  PetscCall(EPSSetFromOptions(cross->eps));
  PetscFunctionReturn(PETSC_SUCCESS);
}
```
1) 
```
PetscOptionsHeadBegin(PetscOptionsObject,"SVD Cross Options");

    PetscCall(PetscOptionsBool("-svd_cross_explicitmatrix","Use cross explicit matrix","SVDCrossSetExplicitMatrix",cross->explicitmatrix,&val,&set));
    if (set) PetscCall(SVDCrossSetExplicitMatrix(svd,val));

  PetscOptionsHeadEnd();
```
Предоставляет пользователю возможность выбрать использование явной матрицы вместо шельковой матрицы через интерфейс опций. Если выбрано использование явной матрицы, то вызывается функция SVDCrossSetExplicitMatrix, которая устанавливает режим работы с явной матрицей.
```
if (!cross->eps) PetscCall(SVDCrossGetEPS(svd,&cross->eps));
  if (!cross->explicitmatrix && !cross->usereps) {
    /* use as default an ST with shell matrix and Jacobi */
    PetscCall(EPSGetST(cross->eps,&st));
    PetscCall(STSetMatMode(st,ST_MATMODE_SHELL));
  }
```
Если пользователь не выбрал явную матрицу и не использовал встроенный EPS, автоматически применяется тип преобразования с использованием режима Shell Matrix. Это означает, что ST наносит ущерб работе с матрицами в режиме настройки, который настраивает для многих задач сингулярного разложения.

# Функция SVDCrossSetExplicitMatrix_Cross(SVD svd,PetscBool explicitmatrix) -изменяет поведение объекта SVDCross таким образом, что он начинает использовать явную матрицу вместо матрици SHELL. Это изменение поведения возможно только в том случае, если текущий режим работы не соответствует ожидаемому значению параметра explicitmatrix
```
static PetscErrorCode SVDCrossSetExplicitMatrix_Cross(SVD svd,PetscBool explicitmatrix)
{
  SVD_CROSS *cross = (SVD_CROSS*)svd->data;

  PetscFunctionBegin;
  if (cross->explicitmatrix != explicitmatrix) {
    cross->explicitmatrix = explicitmatrix;
    svd->state = SVD_STATE_INITIAL;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}
```
```
if (cross->explicitmatrix != explicitmatrix) {
    cross->explicitmatrix = explicitmatrix;
    svd->state = SVD_STATE_INITIAL;
  }
```
Проверяет текущий состояние объекта cross и если флаг explicitmatrix отличается от желаемого значения, то он устанавливается равным желаемому значению, а состояние объекта svd устанавливается в начальное состояние.

# Функция SVDCrossSetExplicitMatrix(SVD svd,PetscBool explicitmat) - Функция проверяет корректность объектов и пытается применить метод изменения режима работы объекта SVD на использование явной матрицы
```
PetscErrorCode SVDCrossSetExplicitMatrix(SVD svd,PetscBool explicitmat)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscValidLogicalCollectiveBool(svd,explicitmat,2);
  PetscTryMethod(svd,"SVDCrossSetExplicitMatrix_C",(SVD,PetscBool),(svd,explicitmat));
  PetscFunctionReturn(PETSC_SUCCESS);
}
```
Класс SVD проверяется на соответствие классу и наличие нужных флагов.
Параметр explicitmat проверяется на наличие и правильность значения.

# Функция SVDCrossGetExplicitMatrix_Cross(SVD svd,PetscBool *explicitmat) - функция флаг explicitmatrix, который определяет использование явной матрицы для решения задачи собственных значений.
```
static PetscErrorCode SVDCrossGetExplicitMatrix_Cross(SVD svd,PetscBool *explicitmat)
{
  SVD_CROSS *cross = (SVD_CROSS*)svd->data;

  PetscFunctionBegin;
  *explicitmat = cross->explicitmatrix;
  PetscFunctionReturn(PETSC_SUCCESS);
}
```
# Функция SVDCrossGetExplicitMatrix(SVD svd,PetscBool *explicitmat) - функция возвращает указатель на явную матрицу, если режим работы объекта SVD установлен на использование явной матрицы. 
```
PetscErrorCode SVDCrossGetExplicitMatrix(SVD svd,PetscBool *explicitmat)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscAssertPointer(explicitmat,2);
  PetscUseMethod(svd,"SVDCrossGetExplicitMatrix_C",(SVD,PetscBool*),(svd,explicitmat));
  PetscFunctionReturn(PETSC_SUCCESS);
}
```
# Функция = PetscErrorCode SVDCrossSetEPS_Cross(SVD svd,EPS eps) - функция реализует замену объекта eps в структуре SVD_Cross на новый объект eps
```
static PetscErrorCode SVDCrossSetEPS_Cross(SVD svd,EPS eps)
{
  SVD_CROSS      *cross = (SVD_CROSS*)svd->data;

  PetscFunctionBegin;
  PetscCall(PetscObjectReference((PetscObject)eps));
  PetscCall(EPSDestroy(&cross->eps));
  cross->eps     = eps;
  cross->usereps = PETSC_TRUE;
  svd->state     = SVD_STATE_INITIAL;
  PetscFunctionReturn(PETSC_SUCCESS);
}
```
# Функция SVDCrossSetEPS(SVD svd,EPS eps) - функция позволяет настроить связь между объектами SVD и EPS для решения задачи собственных значений, обеспечивая правильное взаимодействие между этими объектами.
```
PetscErrorCode SVDCrossSetEPS(SVD svd,EPS eps)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscValidHeaderSpecific(eps,EPS_CLASSID,2);
  PetscCheckSameComm(svd,1,eps,2);
  PetscTryMethod(svd,"SVDCrossSetEPS_C",(SVD,EPS),(svd,eps));
  PetscFunctionReturn(PETSC_SUCCESS);
}
```
Проверка валидовности объектов svd и eps: Проверка класса SVD, идентификатор класса и наличие необходимых флагов. Проверка классов и интерфейсов для объектов eps.
Если проверка прошла успешно, вызывается метод SVDCrossSetEPS_C для настройки связи между объектами SVD и EPS.


# Функция SVDCrossGetEPS_Cross(SVD svd,EPS *eps) - функция реализует настройку метода SVDCross для решения задачи собственных чисел с использованием объекта EPS
```
static PetscErrorCode SVDCrossGetEPS_Cross(SVD svd,EPS *eps)
{
  SVD_CROSS      *cross = (SVD_CROSS*)svd->data;

  PetscFunctionBegin;
  if (!cross->eps) {
    PetscCall(EPSCreate(PetscObjectComm((PetscObject)svd),&cross->eps));
    PetscCall(PetscObjectIncrementTabLevel((PetscObject)cross->eps,(PetscObject)svd,1));
    PetscCall(EPSSetOptionsPrefix(cross->eps,((PetscObject)svd)->prefix));
    PetscCall(EPSAppendOptionsPrefix(cross->eps,"svd_cross_"));
    PetscCall(PetscObjectSetOptions((PetscObject)cross->eps,((PetscObject)svd)->options));
    PetscCall(EPSSetWhichEigenpairs(cross->eps,EPS_LARGEST_REAL));
    PetscCall(EPSMonitorSet(cross->eps,EPSMonitor_Cross,svd,NULL));
  }
  *eps = cross->eps;
  PetscFunctionReturn(PETSC_SUCCESS);
}
```
1)
```
SVD_CROSS      *cross = (SVD_CROSS*)svd->data;
``` 
Cоздается переменная cross типа SVD_CROSS*, которая указывает на данные объекта svd. Тип данных svd преобразуется к типу SVD_CROSS*.
```
 if (!cross->eps) {
```
Cуществует ли объект eps внутри структуры cross
```
PetscCall(EPSCreate(PetscObjectComm((PetscObject)svd),&cross->eps));
```
Создается новый объект типа EPS и сохраняется в структуре cross
```
PetscCall(PetscObjectIncrementTabLevel((PetscObject)cross->eps,(PetscObject)svd,1));
```
Увеличивает количество вкладок, которые PETSCVIEWERASCII выводит для этого объекта, используя tablevel другого объекта. Это должно быть вызвано сразу после создания объекта.
Входные параметры
PetscObjectIncrementTabLevel((PetscObject)cross->ep - любой объект PETSc, в котором мы меняем вкладку
(PetscObject)svd - объект, предоставляющий вкладку, необязательно передайте NULL, чтобы использовать 0 в качестве предыдущего tablevel для obj
1 - приращение, которое добавляется к вкладке старых объектов
```
    PetscCall(EPSSetOptionsPrefix(cross->eps,((PetscObject)svd)->prefix));
```
Функция EPSAppendOptionsPrefix добавляет указанный префикс ("svd_cross_") к текущему списку параметров конфигурации объекта cross->eps.
```
PetscCall(EPSAppendOptionsPrefix(cross->eps,"svd_cross_"));
```
Добавляется дополнительный префикс "svd_cross_" к  cross->eps
```
PetscCall(PetscObjectSetOptions((PetscObject)cross->eps,((PetscObject)svd)->options));
```
Передаются параметры конфигурации от родительского объекта svd новому объекту cross->eps.
```
PetscCall(EPSSetWhichEigenpairs(cross->eps,EPS_LARGEST_REAL));
```
Режим поиска собственных значений для объекта EPS: будут найдены наибольшие вещественные собственные значения.
```
    PetscCall(EPSMonitorSet(cross->eps,EPSMonitor_Cross,svd,NULL));
```
Устанавливается монитор для отслеживания процесса решения задачи собственных значений.

# Функция SVDCrossGetEPS(SVD svd,EPS *eps) - функция, которая помогает получить объект EPS из объекта SVD
```
PetscErrorCode SVDCrossGetEPS(SVD svd,EPS *eps)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscAssertPointer(eps,2);
  PetscUseMethod(svd,"SVDCrossGetEPS_C",(SVD,EPS*),(svd,eps));
  PetscFunctionReturn(PETSC_SUCCESS);
}
```
```
PetscFunctionBegin;
```
```
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
```
Принадлежит ли указанному классу 1 элемент структуры svd 
```
  PetscAssertPointer(eps,2);
```
Не равен ли 2-ой элемент eps, NULL
```
  PetscUseMethod(svd,"SVDCrossGetEPS_C",(SVD,EPS*),(svd,eps));
```
Запрашивает PetscObject для метода, добавленного с помощью PetscObjectComposeFunction(), если он существует, то вызывает его, в противном случае генерирует ошибку.
Входные параметры
svd - объект, например Mat, который не нужно приводить к PetscObject
"SVDCrossGetEPS_C" - имя метода, например, «KSPGMRESSetRestart_C» для функции KSPGMRESSetRestart()
(SVD,EPS*) - типы аргументов для метода, например, (KSP,PetscInt)
(svd,eps) - аргументы для метода, например, (ksp,restart))

# Функция SVDView_Cross(SVD svd,PetscViewer viewer) - функция отвечает за вывод информации о сингулярном разложении (SVD) в формате ASCII через объект PetscViewer
```
static PetscErrorCode SVDView_Cross(SVD svd,PetscViewer viewer)
{
  SVD_CROSS      *cross = (SVD_CROSS*)svd->data;
  PetscBool      isascii;

  PetscFunctionBegin;
  PetscCall(PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii));
  if (isascii) {
    if (!cross->eps) PetscCall(SVDCrossGetEPS(svd,&cross->eps));
    PetscCall(PetscViewerASCIIPrintf(viewer,"  %s matrix\n",cross->explicitmatrix?"explicit":"implicit"));
    PetscCall(PetscViewerASCIIPushTab(viewer));
    PetscCall(EPSView(cross->eps,viewer));
    PetscCall(PetscViewerASCIIPopTab(viewer));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}
```
```
PetscCall(PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii));
```
Функция PetscObjectTypeCompare проверяет, является ли объект viewer типом PETSCVIEWERASCII. Результат сравнения записывается в переменную isascii
```
if (isascii) {
    if (!cross->eps) PetscCall(SVDCrossGetEPS(svd,&cross->eps));
    PetscCall(PetscViewerASCIIPrintf(viewer,"  %s matrix\n",cross->explicitmatrix?"explicit":"implicit"));
    PetscCall(PetscViewerASCIIPushTab(viewer));
    PetscCall(EPSView(cross->eps,viewer));
    PetscCall(PetscViewerASCIIPopTab(viewer));
  }
```
Проверяется наличие объекта EPS в структуре cross. Если его нет, вызывается функция SVDCrossGetEPS для его создания.
Вызывается PetscViewerASCIIPrintf для вывода информации о матрице (явной или неявной).
С помощью PetscViewerASCIIPushTab устанавливается уровень табуляции для вывода.
Вызывается EPSView для отображения содержимого объекта EPS.
С помощью PetscViewerASCIIPopTab восстанавливается предыдущий уровень табуляции.

# Функция SVDReset_Cross(SVD svd) - функция сбрасывает внутреннее состояние объекта SVD, удаляя временные матрицы и подготавливая объект к повторному использованию
```
static PetscErrorCode SVDReset_Cross(SVD svd)
{
  SVD_CROSS      *cross = (SVD_CROSS*)svd->data;

  PetscFunctionBegin;
  PetscCall(EPSReset(cross->eps));
  PetscCall(MatDestroy(&cross->C));
  PetscCall(MatDestroy(&cross->D));
  PetscFunctionReturn(PETSC_SUCCESS);
}
```

```
 PetscCall(EPSReset(cross->eps));
```
Эта функция сбрасывает состояние объекта EPS так, чтобы его можно было использовать заново.
```
  PetscCall(MatDestroy(&cross->C));
  PetscCall(MatDestroy(&cross->D));
```
Удаляются две матрицы, хранимые в полях C и D структуры cross. Эти матрицы могут содержать промежуточные результаты вычислений, которые больше не нужны после сброса.

# Функция SVDDestroy_Cross(SVD svd) - функция удаляет ресурсы, связанные с объектом SVD, включая уничтожение объекта EPS.
```
static PetscErrorCode SVDDestroy_Cross(SVD svd)
{
  SVD_CROSS      *cross = (SVD_CROSS*)svd->data;

  PetscFunctionBegin;
  PetscCall(EPSDestroy(&cross->eps));
  PetscCall(PetscFree(svd->data));
  PetscCall(PetscObjectComposeFunction((PetscObject)svd,"SVDCrossSetEPS_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)svd,"SVDCrossGetEPS_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)svd,"SVDCrossSetExplicitMatrix_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)svd,"SVDCrossGetExplicitMatrix_C",NULL));
  PetscFunctionReturn(PETSC_SUCCESS);
}
```

```
PetscCall(EPSDestroy(&cross->eps));
```
Удаляются данные объекта cross->eps.
```
  PetscCall(PetscFree(svd->data));
```
Освобождается память занимаемая svd->data
```
  PetscCall(PetscObjectComposeFunction((PetscObject)svd,"SVDCrossSetEPS_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)svd,"SVDCrossGetEPS_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)svd,"SVDCrossSetExplicitMatrix_C",NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)svd,"SVDCrossGetExplicitMatrix_C",NULL));
```
Эти строки очищают составные функции, связанные с объектом svd, устанавливая их в NULL.
Так как функция связывает с заданным объектом PETSc - NULL.
PetscObjectComposeFunction(PetscObject obj, const char name[], void (*fptr)(void))
Входные параметры
obj - объект PETSc; он должен быть приведен к (PetscObject), например, PetscObjectCompose((PetscObject)mat,…);
name - имя, связанное с дочерней функцией
fptr - указатель на функцию


# Функция PetscErrorCode SVDCreate_Cross(SVD svd) - функция SVDCreate_Cross создает и инициализирует объект SVD для кросс-сингулярного разложения. Она выделяет память для хранения данных, устанавливает необходимые операции и дополнительные функции для работы с объектом.
```
SLEPC_EXTERN PetscErrorCode SVDCreate_Cross(SVD svd)
{
  SVD_CROSS      *cross;

  PetscFunctionBegin;
  PetscCall(PetscNew(&cross));
  svd->data = (void*)cross;

  svd->ops->solve          = SVDSolve_Cross;
  svd->ops->solveg         = SVDSolve_Cross;
  svd->ops->solveh         = SVDSolve_Cross;
  svd->ops->setup          = SVDSetUp_Cross;
  svd->ops->setfromoptions = SVDSetFromOptions_Cross;
  svd->ops->destroy        = SVDDestroy_Cross;
  svd->ops->reset          = SVDReset_Cross;
  svd->ops->view           = SVDView_Cross;
  svd->ops->computevectors = SVDComputeVectors_Cross;
  PetscCall(PetscObjectComposeFunction((PetscObject)svd,"SVDCrossSetEPS_C",SVDCrossSetEPS_Cross));
  PetscCall(PetscObjectComposeFunction((PetscObject)svd,"SVDCrossGetEPS_C",SVDCrossGetEPS_Cross));
  PetscCall(PetscObjectComposeFunction((PetscObject)svd,"SVDCrossSetExplicitMatrix_C",SVDCrossSetExplicitMatrix_Cross));
  PetscCall(PetscObjectComposeFunction((PetscObject)svd,"SVDCrossGetExplicitMatrix_C",SVDCrossGetExplicitMatrix_Cross));
  PetscFunctionReturn(PETSC_SUCCESS);
}
```
```
svd->ops->solve          = SVDSolve_Cross;
  svd->ops->solveg         = SVDSolve_Cross;
  svd->ops->solveh         = SVDSolve_Cross;
  svd->ops->setup          = SVDSetUp_Cross;
  svd->ops->setfromoptions = SVDSetFromOptions_Cross;
  svd->ops->destroy        = SVDDestroy_Cross;
  svd->ops->reset          = SVDReset_Cross;
  svd->ops->view           = SVDView_Cross;
  svd->ops->computevectors = SVDComputeVectors_Cross;
```
Задаются функции, которые будут выполнять соответствующие операции над объектом SVD. 
```
  PetscCall(PetscObjectComposeFunction((PetscObject)svd,"SVDCrossSetEPS_C",SVDCrossSetEPS_Cross));
  PetscCall(PetscObjectComposeFunction((PetscObject)svd,"SVDCrossGetEPS_C",SVDCrossGetEPS_Cross));
  PetscCall(PetscObjectComposeFunction((PetscObject)svd,"SVDCrossSetExplicitMatrix_C",SVDCrossSetExplicitMatrix_Cross));
  PetscCall(PetscObjectComposeFunction((PetscObject)svd,"SVDCrossGetExplicitMatrix_C",SVDCrossGetExplicitMatrix_Cross));
```
Связывает функцию с заданным объектом PETSc.
Общий вид:
```
PetscObjectComposeFunction(PetscObject obj, const char name[], void (*fptr)(void))
```
Входные параметры
obj - объект PETSc; он должен быть приведен к (PetscObject), например, PetscObjectCompose((PetscObject)mat,…);
name - имя, связанное с дочерней функцией
fptr - указатель на функцию

#4. References.

##### 1. SLEPc Users Manual Scalable Library for Eigenvalue Problem Computations, Carmen Campos Jose E. Roman, Eloy Romero Andres Tomas. Раздел EPS: Eigenvalue Problem Solver, страница 17.
###### 2. http://www.grycap.upv.es/slepc.
###### 3. https://petsc.org/release/
###### 4. https://slepc.upv.es/slepc-main/include/slepc/private/
