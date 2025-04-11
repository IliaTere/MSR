#!/bin/bash

# Цвета для вывода
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
MAGENTA='\033[0;35m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

# Проверка наличия скомпилированных программ
if [ ! -f "residual.out" ] || [ ! -f "error.out" ]; then
    echo -e "${RED}Ошибка: не найдены исполняемые файлы residual.out (метод минимальных невязок) и/или error.out (метод минимальных ошибок).${NC}"
    echo "Пожалуйста, убедитесь, что обе программы скомпилированы и доступны в текущей директории."
    exit 1
fi

# Функция для запуска теста и сравнения результатов
run_comparison() {
    local test_name=$1
    local a=$2
    local b=$3
    local c=$4
    local d=$5
    local nx=$6
    local ny=$7
    local func=$8
    local eps=$9
    local maxit=${10}
    local threads=${11}
    
    echo -e "${MAGENTA}========== ТЕСТ: $test_name ==========${NC}"
    echo -e "${CYAN}Параметры: a=$a b=$b c=$c d=$d nx=$nx ny=$ny func=$func eps=$eps maxit=$maxit threads=$threads${NC}"
    echo
    
    # Запуск программы с методом минимальных невязок
    echo -e "${YELLOW}Метод минимальных невязок:${NC}"
    residual_output=$(./residual.out $a $b $c $d $nx $ny $func $eps $maxit $threads)
    echo "$residual_output"
    
    # Извлечение значений невязок и итераций
    residual_r1=$(echo "$residual_output" | grep -oP "R1 = \K[0-9e\.\-]+")
    residual_r2=$(echo "$residual_output" | grep -oP "R2 = \K[0-9e\.\-]+")
    residual_r3=$(echo "$residual_output" | grep -oP "R3 = \K[0-9e\.\-]+")
    residual_r4=$(echo "$residual_output" | grep -oP "R4 = \K[0-9e\.\-]+")
    residual_it=$(echo "$residual_output" | grep -oP "It = \K[0-9]+")
    residual_t1=$(echo "$residual_output" | grep -oP "T1 = \K[0-9\.]+")
    
    echo
    
    # Запуск программы с методом минимальных ошибок
    echo -e "${GREEN}Метод минимальных ошибок:${NC}"
    error_output=$(./error.out $a $b $c $d $nx $ny $func $eps $maxit $threads)
    echo "$error_output"
    
    # Извлечение значений невязок и итераций
    error_r1=$(echo "$error_output" | grep -oP "R1 = \K[0-9e\.\-]+")
    error_r2=$(echo "$error_output" | grep -oP "R2 = \K[0-9e\.\-]+")
    error_r3=$(echo "$error_output" | grep -oP "R3 = \K[0-9e\.\-]+")
    error_r4=$(echo "$error_output" | grep -oP "R4 = \K[0-9e\.\-]+")
    error_it=$(echo "$error_output" | grep -oP "It = \K[0-9]+")
    error_t1=$(echo "$error_output" | grep -oP "T1 = \K[0-9\.]+")
    
    echo
    
    # Сравнение результатов
    echo -e "${BLUE}Сравнение методов:${NC}"
    echo -e "Невязка R1 (норма невязки): Невязки ${residual_r1} vs Ошибки ${error_r1}"
    echo -e "Невязка R2 (интегральная невязка): Невязки ${residual_r2} vs Ошибки ${error_r2}"
    echo -e "Невязка R3 (максимальная ошибка): Невязки ${residual_r3} vs Ошибки ${error_r3}"
    echo -e "Невязка R4 (интегральная ошибка): Невязки ${residual_r4} vs Ошибки ${error_r4}"
    echo -e "Количество итераций: Невязки ${residual_it} vs Ошибки ${error_it}"
    echo -e "Время выполнения: Невязки ${residual_t1} vs Ошибки ${error_t1}"
    
    # Определение, какой метод лучше
    local wins_residual=0
    local wins_error=0
    
    if (( $(echo "$residual_r1 < $error_r1" | awk '{print ($1 < $3)}') )); then
        wins_residual=$((wins_residual+1))
    else
        wins_error=$((wins_error+1))
    fi
    
    if (( $(echo "$residual_r2 < $error_r2" | awk '{print ($1 < $3)}') )); then
        wins_residual=$((wins_residual+1))
    else
        wins_error=$((wins_error+1))
    fi
    
    if (( $(echo "$residual_r3 < $error_r3" | awk '{print ($1 < $3)}') )); then
        wins_residual=$((wins_residual+1))
    else
        wins_error=$((wins_error+1))
    fi
    
    if (( $(echo "$residual_r4 < $error_r4" | awk '{print ($1 < $3)}') )); then
        wins_residual=$((wins_residual+1))
    else
        wins_error=$((wins_error+1))
    fi
    
    if [ "$residual_it" -lt "$error_it" ]; then
        wins_residual=$((wins_residual+1))
    else
        wins_error=$((wins_error+1))
    fi
    
    if (( $(echo "$residual_t1 < $error_t1" | awk '{print ($1 < $3)}') )); then
        wins_residual=$((wins_residual+1))
    else
        wins_error=$((wins_error+1))
    fi
    
    echo
    echo -e "Количество преимуществ: Метод минимальных невязок ($wins_residual) vs Метод минимальных ошибок ($wins_error)"
    
    if [ "$wins_residual" -gt "$wins_error" ]; then
        echo -e "${YELLOW}Для данного теста метод минимальных невязок показал лучшие результаты.${NC}"
    elif [ "$wins_residual" -lt "$wins_error" ]; then
        echo -e "${GREEN}Для данного теста метод минимальных ошибок показал лучшие результаты.${NC}"
    else
        echo -e "${BLUE}Методы показали одинаковые результаты.${NC}"
    fi
    
    echo -e "${MAGENTA}========== КОНЕЦ ТЕСТА ==========${NC}"
    echo
}

echo "==========================================================="
echo "  СРАВНЕНИЕ МЕТОДА МИНИМАЛЬНЫХ НЕВЯЗОК И МЕТОДА МИНИМАЛЬНЫХ ОШИБОК"
echo "==========================================================="
echo

# Тест 1: Стандартный случай (базовое сравнение)
run_comparison "Стандартный случай" 0 1 0 1 50 50 0 1e-8 1000 4

# Тест 2: Вытянутая область (плохая обусловленность матрицы)
run_comparison "Вытянутая область" -10 10 -0.5 0.5 150 15 3 1e-10 1000 4

# Тест 3: Высокая точность (проверка сходимости)
run_comparison "Высокая точность" 0 1 0 1 100 100 0 1e-12 1000 4

# Тест 4: Сложная функция (проверка работы на сложных функциях)
run_comparison "Сложная функция" 0 1 0 1 80 80 7 1e-8 1000 4

# Тест 5: Большая сетка (проверка скорости и масштабируемости)
run_comparison "Большая сетка" 0 1 0 1 200 200 0 1e-8 1000 8

# Тест 6: Экстремальный случай (проверка стабильности)
run_comparison "Экстремальный случай" -5 5 -0.2 0.2 120 12 7 1e-10 2000 4

echo "==========================================================="
echo "  ЗАВЕРШЕНО СРАВНЕНИЕ МЕТОДОВ"
echo "===========================================================" 