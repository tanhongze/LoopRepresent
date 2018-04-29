Loop Represent系列程序v3
20180228  BY:谭弘泽

文档结构：
--主文件夹
 |--common.h
 |  公共头文件
 |--loop.h
 |  圈表示头文件
 |--worm.h
 |  worm算法头文件
 |--readme.txt
 |  文档结构
 |--算法说明.docx
 |  关于本程序组实现的算法的部分说明
 |   
 |--exhaustion
 | |  枚举法相关程序
 | |--propagator.cpp
 | |    枚举法算传播子预处理
 | |--propagator.py
 | |    枚举法算传播子-绘制物理量
 | |--state.cpp
 | |    枚举法算圈表示
 | |--check.py
 | |    枚举自洽性检查
 | |--plot.py
 |      枚举法算圈表示-相关图像
 |   
 |--graph
 | |  图像相关脚本
 | |--loop
 | |    圈表示-图像相关脚本
 | |--worm
 |      worm算法-python图像相关代码
 |   
 |--initial
 | |  初始化数据文件
 | |--loop_initial.h
 | |    圈表示初始化数据
 | |--worm_initial.h
 |      worm算法初始化数据
 |   
 |--main
 | |  主程序相关文件
 | |--worm_propagators.cpp
 | |  计算worm到达各个位置的相位叠加，进行傅里叶变换，并输出到.mbf文件  
 | |--worm_mc.cpp
 | |  搜索临界质量
 | |--loop_main.cpp
 | |  圈表示主文件
 | |--worm_psibargamma5psi_main.cpp
 | |  输出零动量传播子和拟合后的有效质量
 | |--worm_observers_psibarpsi.cpp
 | |  输出物理量到python数据，可以利用../exhaustion/propagator.py程序的结果进行对比（仅限3*3格点）
 | |--worm_psi_main.cpp
 |    输出单夸克零动量传播子和拟合后的有效质量
 | 
 | 
 |--make_initial
 | |  初始化程序文件
 | |--loop_initial.h
 | |    圈表示初始化数据生成程序
 | |--worm_initial.h
 |      worm算法初始化数据生成程序
 |--other
 | |
 | |--ci.cpp
 | |  计算系数通项公式中的组合计数部分
 | |--context.cpp
 | |  局部信息的一些测试
 | |--windings_sign.cpp
 |    对比文献中求和公式和化简后适用于一般情况的Z0100...00的系数
 | 
 |--utils
 |  工具包，包括实现打印格点信息，计时，傅里叶变换，粗糙的拟合功能，和一些作为常数的类
 |--version
    程序编译运行时自动推进的版本信息
 
 注：编译的时候需要打开-std=C++11选项并链接openmp相关库文件