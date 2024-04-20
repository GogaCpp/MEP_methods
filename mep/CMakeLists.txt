# Создаем библиотеку
add_library(mep mepClass.cpp)

# Указываем заголовочные файлы (опционально, но рекомендуется)
target_include_directories(mep PUBLIC .)
