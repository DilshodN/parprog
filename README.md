# parprog

Библиотека для работы с матрицами с использованием многопоточности от `std::thread`

Код полностью покрыт unit-test-ами (googletest)

**Сборка:**
```
mkdir build
cd build
cmake ../CMakeLists.txt
make
```
**Запуск:**
1. Программы (вывод результатов операции над матрицами при различном числе потоков):
`./build/multithread_matrix`
2. Тестов (googletest)
`./build/test/test`

Результаты приведены в файлике: `results.txt`
**Device:**


OS: macOS Catalina


Processor: 1,1 GHz 2‑x kernel IntelCore m3
