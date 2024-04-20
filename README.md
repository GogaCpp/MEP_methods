# Сравнение численных методов для поиска пути с минимальным перепадом энергии
## Елизаров Андрей Андреевич



При изучении химических и физических процессов, которые включают переход моделируемой системы из одного состояния в другое, возникает задача поиска пути с минимальным перепадом энергии (англ. minimum energy path). В простейшем случае переход между состояниями представляется в виде пути из одного энергетического минимума в другой, который включает в себя потенциальный барьер. Если рассматривать переход при непрерывной поверхности потенциальной энергии, то искомый путь возможно определить в виде кривой,либо множеством точек в n- мерном пространстве

Основная цель моего исследования сравнить два самых популярных численных метода(NEB и string).

В этом репозитории представлена реализация этих двух методов и тестовый пример.



## Запуск тестового примера


```sh
cd build
make
./example
```