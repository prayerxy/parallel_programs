并行计算复习



## 并行概述

假定已有求解问题的串行算法，我们将其改写为并行版本。

### 任务分解

- ==将计算任务进行分解，交由众多进程/线程并发执行==
- 保持依赖关系：计算结果必须与串行算法一致
- 额外开销：锁、同步、barrier等的开销



==竞争条件与数据依赖==：

如果执行结果依赖于多个事件的时序，那么存在竞争条件。

- 数据依赖：内存操作的时序条件，必须保持，以保证结果的正确性
- 同步：将多个线程的执行串行化（比如只有多个线程都完成了某一操作，这些线程才会进行下一步动作）

例子：n个数求和

使用计算任务划分，此时：

- 划分任务，一个线程进行n/t个数求和，最后sum+=此线程的求和值
- 但是由于每个线程在循环内部都会用到sum，而sum是一个全局变量，==任何时候必须只有一个线程对sum操作==
- 思考如何粗粒度加锁/同步
  1. 每一个线程先完成自己的局部求和，最后在循环外将局部和加到全局sum时，对此操作加锁
  2. 消除锁，共享局部数组my_sum[t]，设置同步(每一个线程均完成自己的局部和求解后)，0线程单独进行sum的求和操作
  3. 多核并行求全局和：不是一个线程负责所有线程的局部和相加，而是递归分解

==开销==

- 使用锁mutex、barrier障碍，进行同步
- 划分任务的额外代码，根据每个线程的id划分
- 更多地进行本地运算，少进行全局值的更改
- 可能出现负载不均，每个线程的任务划分要合理



### 数据并行

- 将求解问题涉及到的数据划分给多个核心
- 每个核心对 不同数据进行相似的计算
- 类似于批改试卷，多个人批改试卷，批改的行为相同，而试卷是数据，对数据并行
- ==行/列划分==、==2维划分==



### 其他任务划分

#### 递归分解

- 分治策略
- 问题求解，转化为子问题的求解，子问题求解继续分而治之，分解停止——原子问题。
- 但是枢轴要选取恰当，不然容易负载不均



#### 搜索分解

与数据分解区别：一旦一个任务找到解，全部任务可以停止

#### 通信优化

- 更多地访问局部数据/最近访问过的数据，类似于cache
- 最小化数据交互量：中间结果保存为局部数据，最终结果计算时才通信
- 最小化冲突：重构算法，避免多个任务同时访问共享数据，多个进行同时向一个进程发送消息等
- 最小化交互频率：重构算法，大块访问数据

#### 矩阵相乘示例

C=AB

划分子矩阵，每个进程负责C的一个子矩阵运算，注意分配任务与进程标号结合，使得每次子矩阵相乘用到的子矩阵不同，错开。

<img src="C:/Users/%E5%8D%8E%E7%9B%96%E5%B0%86%E5%80%BE/AppData/Roaming/Typora/typora-user-images/image-20240608104536636.png" alt="image-20240608104536636" style="zoom:67%;" />

#### 通信与通信重叠

<img src="C:/Users/%E5%8D%8E%E7%9B%96%E5%B0%86%E5%80%BE/AppData/Roaming/Typora/typora-user-images/image-20240608104616199.png" alt="image-20240608104616199" style="zoom:67%;" />

传递四个消息，可以把a做四次；但是不如c，c类似于流水线，每个进程传播多个消息，通信重叠度大



### 分析并行算法

1. 串行算法评价：时间复杂度

2. 并行算法评价：输入规模、处理器数目、通信速度

3. 评价标准：运行事件；加速比

4. 并行程序设计的复杂性：

   1. 足够的并发度(Amdahl定律)
   2. 并发粒度：独立的计算任务的大小
   3. 局部性：对临近数据计算
   4. 负载均衡
   5. 协调和同步

5. 并行算法额外开销

   1. 进程间通信：最大开销

   2. 进程空闲：负载不均

   3. 额外计算：需要同步

   4. $T_o$=$pT_p -T_s$，其中$T_P$是并行算法时间，p是核的数量

   5. 加速比$S=T_s/T_p$，一般S小于等于p，但是S可能大于p，超线性：串行算法计算量>并行算法；硬件问题不利于串行算法

      - 数据量较大，无法放入cache
      - 但是并行后数据量变小，每个处理器的数据量可以放入cache，性能提高

   6. 代价最优：Ts与pTp同阶

      - 如果不是代价最优，当p固定，问题规模n变大时，效率会越来越低
      - E=Ts/pTp    cost=pTp    E<=1

   7. Adamhl定律：随着p增大（核数量变大），S趋于饱和，E减小(开销变大)

      - p不变，n增大，S和E都会增大。

   8. 等效率函数  $W=T_s$，评估并行运行时间

      $T_P = \frac{W + T_o(W, p)}{p}$

      $$S = \frac{W}{T_P} = \frac{Wp}{W + T_o(W, p)}$$

      $E = \frac{S}{p} = \frac{W}{W + T_o(W, p)} = \frac{1}{1 + \frac{T_o(W, p)}{W}}$   其中To代表并行算法总额外开销



### 并行算法实现 SYCL/DPC++

异构编程

SYCL（Single-source C++ for Heterogeneous Computing）是一个开放标准，由Khronos Group定义，允许开发者使用现代C++在异构计算设备上编写并行代码。SYCL 的特点是单一源代码，即主机代码和设备代码可以写在同一个文件中。

==DPC++==（Data Parallel C++）则是 Intel 在 SYCL 基础上的扩展，作为 ==oneAPI==的核心编程模型。

<img src="C:/Users/%E5%8D%8E%E7%9B%96%E5%B0%86%E5%80%BE/AppData/Roaming/Typora/typora-user-images/image-20240608112212607.png" alt="image-20240608112212607" style="zoom: 80%;" />

#### 1.SYCL缓冲方法

使用 `buffer` 和 `accessor`

`buffer` 和 `accessor` 是 SYCL 中管理和传输数据的一种常用方法。`buffer` 负责在主机和设备之间传输数据，而 `accessor` 允许设备内核访问缓冲区中的数据。在这种方式下，SYCL 隐式管理数据同步和传输。

```c++
#include <CL/sycl.hpp>
using namespace sycl;
constexpr int N = 16;
int main() {
    std::vector<double> v(N, 10);
    queue q;
    // 定义一个作用域，用于管理 buffer 的生命周期
    {
        buffer buf(v);
        q.submit([&](handler &h) {
            accessor a(buf, h);
            h.parallel_for(N, [=](auto i) {
                a[i] -= 2;
            });
        });
        // buffer buf 在这里会被销毁，隐式同步
    }

    // 现在数据已经同步回主机，可以安全地访问
    for (int i = 0; i < N; i++) {
        std::cout << v[i] << "\n";
    }
    return 0;
}

```

创建一个缓冲区 `buffer buf(v);` 时，SYCL 会将向量 `v` 的数据复制到缓冲区 `buf` 中，并管理该数据在设备和主机之间的同步。而当缓冲区被销毁时，SYCL 运行时会确保所有与缓冲区相关的任务都已完成。这包括将设备上的数据复制回主机，使主机上的向量 `v` 包含最新的计算结果。

```c++
int main() {
    // Host memory buffer that device will write data back before destruction.
    float(*c_back)[P] = new float[M][P];
    
    queue q(default_selector_v);
    buffer<float, 2> a_buf(range(M, N));
    buffer<float, 2> b_buf(range(N, P));
    buffer c_buf(reinterpret_cast<float *>(c_back), range(M, P));
    //该缓冲区与主机内存中的c_back绑定。
    clock_t start = clock();
    
    q.submit([&](auto &h) {
        accessor a(a_buf, h, write_only);
        h.parallel_for(range(M, N), [=](auto index) { a[index] = 1.0f; });
    });
    
    q.submit([&](auto &h) {
        accessor b(b_buf, h, write_only);
        h.parallel_for(range(N, P), [=](auto index) { b[index] = 1.0f; });
    });
    
    q.submit([&](auto &h) {
        accessor a(a_buf, h, read_only);
        accessor b(b_buf, h, read_only);
        accessor c(c_buf, h, write_only);
        int width_a = a_buf.get_range()[1];
        h.parallel_for(range(M, P), [=](auto index) {
            int row = index[0];
            int col = index[1];
            float sum = 0.0f;
            for (int i = 0; i < width_a; i++) {
                sum += a[row][i] * b[i][col];
            }
            c[index] = sum;
        });
    });
    
    clock_t finish = clock();
    cout << "Run time:" << (double)(finish - start) / CLOCKS_PER_SEC << "\n";
}

```

#### 2.基于指针的方法USM

隐式内存传输`malloc_shared` 

`malloc_shared` 是另一种方式，通过分配共享内存，可以在主机和设备之间直接共享数据，而无需显式的数据传输。`malloc_shared` 分配的内存区域在主机和设备上都可见，且内存同步由 SYCL 运行时自动处理。

USM基于指针的替代，而无需将所有数据变成buffer缓冲区。

```c++
#include <CL/sycl.hpp>
using namespace sycl;

constexpr int N = 16;

int main() {
    queue q;
    double* data = malloc_shared<double>(N, q);

    // 初始化数据
    for (int i = 0; i < N; ++i) {
        data[i] = 10;
    }

    // 提交并行任务
    q.parallel_for(N, [=](auto i) {
        data[i] -= 2;
    }).wait();  //wait同步

    // 主机端访问数据
    for (int i = 0; i < N; ++i) {
        std::cout << data[i] << "\n";
    }

    // 释放共享内存
    free(data, q);

    return 0;
}
```



## SIMD

**概念：多个算术运算->一个SIMD操作**    多个取数/存结果操作->一个更宽的内存操作

主要是针对循环内部的算术

- x86平台 SSE/AVX256/AVX512
- ARM平台 Neon
- 向量化
  - 标量： add r1,r2,r3
  - simd：vadd<sws> v1,v2,v3
- 指令集扩展：编程接口类似于函数调用

### SIMD编程的复杂性

- 底层编程：使用intrinsic；高层编程：使用编译器，而不是自己显式调用
- 数据必须对齐，在内存中要连续存储
- 打包/解包开销：
  - 打包源运算对象——拷贝至连续内存区域
  - 解包目的运算对象——拷贝回内存
- 对齐的内存访问
  - 地址总是向量长度的倍数（例如16字节）
  - 只有一个超级字的读/写操作
- 未对齐的内存访问
  - 地址不是16字节的整数倍
  - 有时硬件会帮做
  - 可能会出现性能下降，甚至结果错误



### x86平台SSE与AVX

1. SSE128位，支持16个8位整数同时运算、8个16位整数
2. AVX256，扩展SSE至256位
3. AVX512，扩展至512位
4. 对应于C/C++的intrinsics，编译器能识别的函数

```c++
//数据移动和相关初始化
_m128 va=_mm_load_ps(addr)   //从此位置取n位至寄存器中
//算术运算
__mm_div_ps(va,vb) //等

//最后从向量寄存器中拷贝回内存
_mm_store_ps(m[k] + j, va)
```

注意load和loadu不一样，load是严格对齐，loadu可以不对齐

#### 高斯消去simd

```c++
for (int k = 0; k < N; k++)
    {
        for (int j = k + 1; j < N; j++)
            m[k][j] = m[k][j] / m[k][k];
        m[k][k] = 1.0;
        for (int i = k + 1; i < N; i++)
        {
            for (int j = k + 1; j < N; j++)
                m[i][j] = m[i][j] - m[i][k] * m[k][j];
            m[i][k] = 0;
        }
    }
```

注意下面的对齐版本的AVX，首先检测地址是否32字节对齐(一个元素是4字节，检测下标是否为8的倍数即可)

```c++
for (int k = 0; k < N; k++) {
    __m256 vt = _mm256_set1_ps(m[k][k]);//用m[k][k]初始化向量vt
    int j = k + 1;
    if (j % 8 != 0) {
        while (j % 8 != 0) {
            m[k][j] = m[k][j] / m[k][k];
            j++;
        }
    }
    for (j; j + 8 <= N; j += 8) {
        __m256 va = _mm256_load_ps(m[k] + j);
        va = _mm256_div_ps(va, vt);
        _mm256_store_ps(m[k] + j, va);//从向量寄存器储存到内存
    }
    for (j; j < N; j++) {
        m[k][j] = m[k][j] / m[k][k];
    }
    m[k][k] = 1.0;
    for (int i = k + 1; i < N; i++) {
        __m256 vaik = _mm256_set1_ps(m[i][k]);
        j = k + 1;
        if (j % 8 != 0) {
            while (j % 8 != 0) {
                m[i][j] = m[i][j] - m[i][k] * m[k][j];
                j++;
            }
        }
        for (j; j + 8 <= N; j += 8) {//内层展开，对j2
            __m256 vakj = _mm256_load_ps(m[k] + j);
            __m256 vaij = _mm256_load_ps(m[i] + j);
            vakj = _mm256_mul_ps(vaik, vakj);
            vaij = _mm256_sub_ps(vaij, vakj);
            _mm256_store_ps(m[i] + j, vaij);
        }
        for (j; j < N; j++) {
            m[i][j] = m[i][j] - m[i][k] * m[k][j];
        }
        m[i][k] = 0;
    }
}
```



### arm平台Neon

- 128位 4个32位int型整数运算

- 常用指令：VLD、VST

- ```c++
  vld1q_f32   //1为非交错模式
  vmulq_f32
  vsubq_f32
  vst1q_f32
  ```

  

#### 高斯消去simd

```c++
 for (int k = 0; k < N; k++) {
    float32x4_t vt = vmovq_n_f32(m[k][k]);//用m[k][k]初始化向量vt
    int j = k + 1;
    if (j % 4 != 0) {
        while (j % 4 != 0) {
            m[k][j] = m[k][j] / m[k][k];
            j++;
        }
    }
    for (j; j + 4 <= N; j += 4) {
        float32x4_t va = vld1q_f32(m[k] + j);
        va = vdivq_f32(va, vt);
        vst1q_f32(m[k] + j, va);//从向量寄存器储存到内存
    }
    for (j; j < N; j++) {
        m[k][j] = m[k][j] / m[k][k];
    }
    m[k][k] = 1.0;
    for (int i = k + 1; i < N; i++) {
        float32x4_t vaik = vmovq_n_f32(m[i][k]);
        j = k + 1;
        if (j % 4 != 0) {
            while (j % 4 != 0) {
                m[i][j] = m[i][j] - m[i][k] * m[k][j];
                j++;
            }
        }
        for (j; j + 4 <= N; j += 4) {//内层展开，对j2
            float32x4_t vakj = vld1q_f32(m[k] + j);
            float32x4_t vaij = vld1q_f32(m[i] + j);
            vakj = vmulq_f32(vaik, vakj);
            vaij = vsubq_f32(vaij, vakj);
            vst1q_f32(m[i] + j, vaij);
        }
        for (j; j < N; j++) {
            m[i][j] = m[i][j] - m[i][k] * m[k][j];
        }
        m[i][k] = 0;
    }
}
```



### 特殊高斯消去算法

外层循环是每一个消元子，内层循环是消元行；如果当前消元子不为空，对所有消元行消元即可；如果当前消元子为空，从被消元行0下标开始查找，是否有可以升格的被消元行，如果有，那么升格，再继续对余下的所有被消元行消元。



```c++
void elimination() {
    for (int i = RN - 1; i >= 0; i--) { //消元子
        if (!eliminer_ifnull(i)) {
            for (int j = 0; j < E_LineN; j++) { //对第j个被消元行进行消元
                if (eline[j].num == i) {
                    int ss = 0;
                    for (ss; ss <= i / 8; ss++)  //第i/8个可能对不需要操作的位进行了异或
                        eline[j].bit[ss] ^= eliminer[i][ss];
                    if (eline_ifnull(j)) {
                        eline[j].ifUprade = true;
                        eline[j].num = -2;
                    }
                    else {//重置num
                        for (int n = eline[j].num; n >= 0; n--) {//找到首个1
                            if (eline[j].bit[n / 8] & (1 << (n % 8))) {  //n%8为0-7的位置，n/8为char的位置
                                eline[j].num = n;
                                break;
                            }
                        }
                    }
                }

            }
        }
        else {


            for (int j = 0; j < E_LineN; j++) {
                if (eline[j].num == i) {
                    for (int ss = eline[j].num / 8; ss >= 0; ss--)  //赋值过来
                        eliminer[eline[j].num][ss] = eline[j].bit[ss];
                    eline[j].ifUprade = true;
                    flagUpgrade = 1;
                    nextstart = j;//之前的被消元行的首位不等于i,不用消元
                    eline[j].num = -2;
                    break;
                }
            }

            if (flagUpgrade == 0)
                continue;
            int newid = nextstart + 1;
            for (int j = newid; j < E_LineN; j++) {
                if (eline[j].num == i) {
                    int ss = 0;
                    for (ss; ss <= i / 8; ss++)  //第i/8个可能对不需要操作的位进行了异或
                        eline[j].bit[ss] ^= eliminer[i][ss];
                    if (eline_ifnull(j)) {
                        eline[j].ifUprade = true;
                        eline[j].num = -2;
                    }
                    else {//重置num
                        for (int n = eline[j].num; n >= 0; n--) {//找到首个1
                            if (eline[j].bit[n / 8] & (1 << (n % 8))) {  //n%8为0-7的位置，n/8为char的位置
                                eline[j].num = n;
                                break;
                            }
                        }
                    }
                }
            }
        }

    }
}
```

simd很简单，直接对所有异或操作进行并行即可

```c++
 for (ss; ss + 64 <= i / 8; ss += 64) { //32个字节，8个int
    __m512i ve = _mm512_loadu_epi8(eline[j].bit + ss);
    __m512i vr = _mm512_loadu_epi8(eliminer[i] + ss);
    __m512i tmp = _mm512_xor_si512(ve, vr);
    _mm512_storeu_epi8(eline[j].bit + ss, tmp);
}

```

- load至向量寄存器
- 然后向量运算异或
- 最后储存至原数据即可

## Pthread

多线程Pthread编程

- 动态线程
  - 主线程等待计算工作，fork新线程分配工作，工作线程完成任务后结束
  - ==创建线程，销毁线程==非常耗时
- 静态线程
  - 创建线程池，线程数量不变，直至整个程序结束
  - 性能更优，但是会浪费资源

线程数据共享

- 全局变量都是共享的

- 函数参数线程id，创建一个线程数据结构

- 线程函数内部定义的变量是私有的

  ```c++
      //初始化barrier
      pthread_barrier_init(&barrier_eliminerNULL, NULL, NUM_THREADS);
      pthread_barrier_init(&barrier_AfterEliminerNull, NULL, NUM_THREADS);
  
      //创建线程
      pthread_t* handles = new pthread_t[NUM_THREADS];
      threadParam_t* param = new threadParam_t[NUM_THREADS];
  
      //创建工作线程
      for (int t_id = 0; t_id < NUM_THREADS; t_id++) {
          param[t_id].t_id = t_id;
          pthread_create(&handles[t_id], NULL, threadFunc, (&param[t_id]));
      }
  
      //销毁进程
      for (int t_id = 0; t_id < NUM_THREADS; t_id++)
          pthread_join(handles[t_id], NULL);
  
  
  	 void pthread_exit(void*value_ptr);//将value_ptr返回给调用者
  ```

### 同步机制

1. 忙等待  while(){}

2. 显式同步：互斥量(锁)   可能导致死锁问题

   ```
   pthread_mutex_lock(&mutex)
   ....
   pthread_mutex_unlock(&mutex)
   ```

   

3. 信号量：sem_wait(sem_t*sem)   sem_post()   初始为0，wait一次信号量-1，post对信号量+1，大于等于0时恢复

4. 使用barrier同步（最简单）

   ```c++
   pthread_barrier_init(&barrier_Divsion, NULL, NUM_THREADS);
   pthread_barrier_init(&barrier_Elimination, NULL, NUM_THREADS);
   
    //第一个同步点
   pthread_barrier_wait(&barrier_Divsion);
   //第二个同步点
   pthread_barrier_wait(&barrier_Elimination);  
   ```

5. mutex底层需要XCHG指令，实现Test-and-Set操作

## OpenMP

- 是pthread的常见替代，更简单
- 可移植与扩展
- 依赖于编译指示，编译器做相关并行化

openmp执行模型：

- ==Fork-join并行执行模型==
- 执行一开始是单进程（主线程）
- 并行结构开始：主线程创建一组线程（工作线程）
- 并行结构结束时：隐式Barrier，线程组同步
- 之后只有主线程继续执行



- Pthread
  - 全局作用域变量共享
  - 栈中分配变量私有
- Openmp
  - 默认shared
  - 可自己指定，shared变量是共享的；private变量是私有的



临界区指令：

- ```
  #pragma omp critical
  ```



规约操作：

- openmp支持规约，避免自己写omp critical等

- ```c++
  sum=0;
  #pragma omp parallel for reduction(+:sum)
  for(i=0;i<100;i++){
  	sum+=array[i];
  }
  ```

  

### 局限

- 并非所有元素级的循环都可以并行化

  #pragma omp parallel for

  - 循环不变量：带符号整数
  - 中止检测，必须是<=,<,>,>=等
  - 每步迭代必须是一个循环不变量
  - 循环体必须无控制流

- ==避免线程创建与销毁开销==，先在外面创建线程，内部再使用omp for并行

  ```c++
  
      int i, j, ss, n, newid;
  #pragma omp parallel num_threads(thread_nums),private(i,j,ss,n,newid)
      {
          for (i = RN - 1; i >= 0; i--) { //消元子
              if (!eliminer_ifnull(i)) {
  #pragma omp for 
                  for (j = 0; j < E_LineN; j++) { //对第j个被消元行进行消元
                      if (eline[j].num == i) {
                          ss = 0;
                          for (ss; ss <= i / 8; ss++)  //第i/8个可能对不需要操作的位进行了异或
                              eline[j].bit[ss] ^= eliminer[i][ss];
                          if (eline_ifnull(j)) {
                              eline[j].ifUprade = true;
                              eline[j].num = -2;
                          }
                          else {//重置num
                              for (n = eline[j].num; n >= 0; n--) {//找到首个1
                                  if (eline[j].bit[n / 8] & (1 << (n % 8))) {  //n%8为0-7的位置，n/8为char的位置
                                      eline[j].num = n;
                                      break;
                                  }
                              }
                          }
                      }
  
                  }
              }
  ```

### 循环分配

schedule子句确定如何在线程间划分循环

- static([chunk]) 分配给每个线程各chunk步迭代；所有线程分配完之后，如果迭代不没有分配完，继续上面的操作
- dynamic([chunk]) 分配给每个线程chunk步迭代，一个线程完成任务之后，再为其分配[chunk]步迭代
  - 动态划分效果可能会更好，达到==负载均衡==，先完成的线程先得到余下的任务
- guided([chunk])  [chunk]呈指数减小
- 默认是static，同时chunk=ceil(iterations/threads_num)



### 一些策略

- 局部性：
  - 减少通信
  - 重用与局部性
    - 数据重用：
      - 多次使用相同或者邻近数据
    - 数据局部性:
      - 重排循环顺序，比如C++数组行主顺序保存
  - 循环展开

### 新特性

- omp single
- 任务并行：同时执行不同计算任务
  - 使用omp sections，下面定义一些omp section区域
  - <img src="C:/Users/%E5%8D%8E%E7%9B%96%E5%B0%86%E5%80%BE/AppData/Roaming/Typora/typora-user-images/image-20240608201024796.png" alt="image-20240608201024796" style="zoom:67%;" />
- omp simd 自动simd向量化

## MPI

MPI 通常假定数据是连续的内存块，因此传输不连续的内存块会导致数据不一致或传输错误，所以在实验时数组不能是动态new出来的



mpi是多进程，由于进程之间地址空间独立，必须要传递消息来通信；而多线程的通信代价显然更小

- 进程之间的通信包括同步与数据移动
- 进程之间不会有数据竞争，但是通信同步可能出错

```c++
MPI_init()
....
MPI_Finalize()
    
    
MPI_Comm_rank(MPI_COMM_WORLD, &myid); //获取id
MPI_Comm_size(MPI_COMM_WORLD, &numprocs);//获取总进程数
```

- 进程与通信域 MPI_COMM_WORLD

- **MPI_Send 阻塞性**：`MPI_Send` 会等待消息数据被复制到系统缓冲区或目标进程开始接收，然后继续执行接下来的代码。这并不意味着目标进程已经处理了该消息。

  ```c++
  int MPI_Send(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm);
  
  ```

  **MPI_Recv 阻塞性**：`MPI_Recv` 会一直等待，直到消息数据完全接收，然后继续执行接下来的代码。

  ```c++
  int MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status *status);
  ```

  

- MPI_ANY_SOURCE、MPI_ANY_TAG   不作筛选

- 可以使用Isend，IRecv等不阻塞的发送与接收

#### MPI编程模型

- 阻塞通信模式：
  - 可能会造成死锁，某个进程发送消息、接收消息
  - Send-receive实际上完成了数据传输与同步操作
  - 发送与接收错开，避免死锁 ：即一个进程是发送、接收，那么另一个进程就是接收、发送。
- 消息传递模型：
  - 对等式
  - 主从式
- 组通信、广播与规约
  - MPI_Bcast，root对于count个进程广播





## GPU

GPU擅长对于大量元素进行相同计算、众核，高并发



### CUDA编程模型

- host：CPU端的代码
- device：GPU端的代码
- Kernel：从host调用，在device端运行的函数

<img src="C:/Users/%E5%8D%8E%E7%9B%96%E5%B0%86%E5%80%BE/AppData/Roaming/Typora/typora-user-images/image-20240608204351786.png" alt="image-20240608204351786" style="zoom:67%;" />

- Thread工作模式
  - 所有线程执行相同的Kernel代码
  - 每个线程使用自己的编号计算不同的数据（内存地址）及执行不同的分支
  - 线程块：Block内线程协同计算，共享内存，原子操作，同步机制；但是Block之间不能协作
  - Block ids，Thread ids。

![image-20240608204628697](C:/Users/%E5%8D%8E%E7%9B%96%E5%B0%86%E5%80%BE/AppData/Roaming/Typora/typora-user-images/image-20240608204628697.png)

- Grid是一组线程块，共享全局内存中的数据，运行时动态调度

例子：==向量相加程序vecAdd的host端==

- 首先在设备商分配A,B,C的空间；将Host中的A,B拷贝至设备上
- 然后Kernel launch code，在设备上运行
- 最后从设备内存中将C拷贝回Host

CUDA的内存层次介绍：

1. Device代码可以：
   - 读/写每线程独占的寄存器
   - 读/写每Kernel共享的全局内存
2. Host代码可以：
   - 在CPU主存与Device全局内存间传输数据
   - <img src="C:/Users/%E5%8D%8E%E7%9B%96%E5%B0%86%E5%80%BE/AppData/Roaming/Typora/typora-user-images/image-20240608205054687.png" alt="image-20240608205054687" style="zoom:67%;" />
   - 数据传输API：
     - cudaMemcpy：Host-Device内存数据传输
     - cudaMalloc在设备端创建数据缓冲区
     - 也可以直接使用`cudaMallocManaged`函数用于在设备和主机之间分配统一内存，分配的内存可以被设备和主机共享

GPU编程例子

```c++
#include <iostream>

// Kernel 函数，用于在 GPU 上执行矢量相加
__global__ void vectorAdd(const float *A, const float *B, float *C, int numElements) {
    //一个block的线程数量*block索引+block内部thread的索引=全局thread的索引
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i < numElements) {
        C[i] = A[i] + B[i];
    }
}

int main() {
    int numElements = 10000; // 矢量元素数量

    // 分配主机上的内存
    float *h_A = new float[numElements];
    float *h_B = new float[numElements];
    float *h_C = new float[numElements];

    // 初始化 A 和 B
    for (int i = 0; i < numElements; ++i) {
        h_A[i] = rand() / (float)RAND_MAX;
        h_B[i] = rand() / (float)RAND_MAX;
    }

    // 分配设备上的内存
    float *d_A, *d_B, *d_C;
    cudaMalloc(&d_A, numElements * sizeof(float));
    cudaMalloc(&d_B, numElements * sizeof(float));
    cudaMalloc(&d_C, numElements * sizeof(float));

    // 将数据从主机复制到设备
    cudaMemcpy(d_A, h_A, numElements * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_B, h_B, numElements * sizeof(float), cudaMemcpyHostToDevice);

    // 启动核函数
    int blockSize = 256;
    int gridSize = (numElements + blockSize - 1) / blockSize;
    //<<<网格中包含块的数量，一个块的线程数量>>>
    vectorAdd<<<gridSize, blockSize>>>(d_A, d_B, d_C, numElements);

    // 将结果从设备复制回主机
    cudaMemcpy(h_C, d_C, numElements * sizeof(float), cudaMemcpyDeviceToHost);

    // 验证结果
    for (int i = 0; i < numElements; ++i) {
        if (fabs(h_A[i] + h_B[i] - h_C[i]) > 1e-5) {
            std::cerr << "Error: Element " << i << " does not match!" << std::endl;
            break;
        }
    }

    // 释放内存
    delete[] h_A;
    delete[] h_B;
    delete[] h_C;
    cudaFree(d_A);
    cudaFree(d_B);
    cudaFree(d_C);

    return 0;
}
```

CUDA变量声明：

- \_\_device\_\_
  - 位于全局内存
  - 与程序一样的生命期
  - Grid中所有的线程均可访问
  - Host通过运行时库可以访问
- \_\_share\_\_
  - 位于某个线程块的共享内存
  - 只有Block内部的所有线程可以访问
- \_\_global\_\_定义一个Kernel函数
- \_\_constant\_\_
  - 位于常量内存
  - 与程序一样的生命期
  - Grid中所有线程均可访问
  - Host通过运行时库可以访问



### CUDA优化

- 全局内存最慢：输入与最终输出数据
- 充分利用共享内存
  - 将数据划分为可放入共享内存的子集
  - 一个Block处理一个子集
    - 数据子集全局内存->共享内存，多线程并发传输
    - 计算结果共享内存->全局内存
- 同步API
  - __syncthreads()  同一个块内的线程同步