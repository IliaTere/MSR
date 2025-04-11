#!/bin/bash

# Цвета для вывода
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Базовый порог невязки
BASE_THRESHOLD=1e-5

# Функция для запуска теста с заданными параметрами
run_test() {
    local a=$1
    local b=$2
    local c=$3
    local d=$4
    local nx=$5
    local ny=$6
    local func=$7
    local eps=$8
    local maxit=$9
    local threads=${10}
    
    # Рассчитываем порог в зависимости от заданной точности
    # Если eps меньше 1e-8, используем порог 1e-5, иначе умножаем eps на 1000
    if (( $(echo "$eps < 1e-8" | awk '{print ($1 < $3)}') )); then
        THRESHOLD=$BASE_THRESHOLD
    else
        THRESHOLD=$(echo "$eps * 1000" | awk '{printf "%.2e", $1 * $3}')
    fi
    
    echo -e "${YELLOW}Тест: a=$a b=$b c=$c d=$d nx=$nx ny=$ny func=$func eps=$eps maxit=$maxit threads=$threads${NC}"
    echo "Порог невязки для этого теста: $THRESHOLD"
    
    # Запуск программы и сохранение вывода
    output=$(./a.out $a $b $c $d $nx $ny $func $eps $maxit $threads)
    
    # Извлечение значений невязок с помощью регулярных выражений
    r1=$(echo "$output" | grep -oP "R1 = \K[0-9e\.\-]+")
    r2=$(echo "$output" | grep -oP "R2 = \K[0-9e\.\-]+")
    r3=$(echo "$output" | grep -oP "R3 = \K[0-9e\.\-]+")
    r4=$(echo "$output" | grep -oP "R4 = \K[0-9e\.\-]+")
    iterations=$(echo "$output" | grep -oP "It = \K[0-9]+")
    
    # Проверка невязок без использования bc
    # Преобразуем научную нотацию в число и сравниваем с порогом
    if [[ "$(echo "$r1 < $THRESHOLD" | awk '{print ($1 < $3)}')" == "1" && \
          "$(echo "$r2 < $THRESHOLD" | awk '{print ($1 < $3)}')" == "1" && \
          "$(echo "$r3 < $THRESHOLD" | awk '{print ($1 < $3)}')" == "1" && \
          "$(echo "$r4 < $THRESHOLD" | awk '{print ($1 < $3)}')" == "1" ]]; then
        echo -e "${GREEN}ТЕСТ ПРОЙДЕН${NC}"
        echo "Невязки: R1=$r1 R2=$r2 R3=$r3 R4=$r4"
        echo "Количество итераций: $iterations"
        return 0
    else
        echo -e "${RED}ТЕСТ НЕ ПРОЙДЕН${NC}"
        echo "Невязки: R1=$r1 R2=$r2 R3=$r3 R4=$r4"
        echo "Количество итераций: $iterations"
        echo "Порог: $THRESHOLD"
        echo "Одна или несколько невязок превышают порог"
        return 1
    fi
}

# Проверка наличия исполняемого файла
if [ ! -f a.out ]; then
    echo -e "${RED}Ошибка: файл a.out не найден. Пожалуйста, скомпилируйте программу.${NC}"
    exit 1
fi

echo "========== НАЧАЛО ТЕСТИРОВАНИЯ =========="
echo "Базовый порог невязки: $BASE_THRESHOLD (для высокой точности eps < 1e-8)"
echo "Для меньшей точности используется порог = eps * 1000"
echo

# Счетчики пройденных и непройденных тестов
passed=0
failed=0

# Тесты с разными функциями и одним потоком
for func in 0 1 2 3; do
    if run_test 0 1 0 1 20 20 $func 1e-8 1000 1; then
        ((passed++))
    else
        ((failed++))
    fi
    echo
done

# Тесты с разными размерами сетки
for size in 10 30 50; do
    if run_test 0 1 0 1 $size $size 0 1e-8 1000 1; then
        ((passed++))
    else
        ((failed++))
    fi
    echo
done

# Тесты с разной точностью
for eps in 1e-6 1e-8 1e-10; do
    if run_test 0 1 0 1 20 20 0 $eps 1000 1; then
        ((passed++))
    else
        ((failed++))
    fi
    echo
done

# Тесты с разным числом потоков
for threads in 1 2 4; do
    if run_test 0 1 0 1 30 30 0 1e-8 1000 $threads; then
        ((passed++))
    else
        ((failed++))
    fi
    echo
done

# Тест с прямоугольной областью
if run_test -1 1 -2 2 20 30 0 1e-8 1000 1; then
    ((passed++))
else
    ((failed++))
fi
echo

echo "========== ИТОГИ ТЕСТИРОВАНИЯ =========="
echo -e "Всего тестов: $((passed + failed))"
echo -e "${GREEN}Пройдено: $passed${NC}"
echo -e "${RED}Не пройдено: $failed${NC}"

if [ $failed -eq 0 ]; then
    echo -e "${GREEN}Все тесты успешно пройдены! Программа работает корректно.${NC}"
    exit 0
else
    echo -e "${RED}Некоторые тесты не пройдены. Необходимо проверить корректность работы программы.${NC}"
    exit 1
fi 