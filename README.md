NADYA project

# GammaCube

## Command Line Interface (CLI)

Проект поддерживает запуск как в интерактивном режиме с визуализацией, так и в батчевом режиме через macro-файл.

### Основные параметры

- `-i, --input`  
  Путь к входному **macro-файлу** с командами для батчевого режима (например, `run.mac`).  
  Если параметр **не указан**, приложение запускается в режиме визуализации.  
  Если параметр **указан**, флаги `-noUI`, `-vd` и `--view-deg` игнорируются.

- `-o, --output-file`  
  Имя выходного `.root` файла.  
  По умолчанию: `GammaCube`.

- `-t, --threads`  
  Количество потоков для многопоточного режима.  
  По умолчанию используется максимальное доступное число ядер.

- `--bins`  
  Количество бинов в выходных гистограммах.

- `--save-seeds`  
  Сохранение зерен в папке event.
  ./Nika -i ../run.mac -fd vertical_down --save-seeds
  ./Nika --replay-event 690

### Параметры сцинтилляции и порогов

- `-ys, --yield-scale`  
  Масштабирование параметра `SCINTILLATIONYIELD` для сцинтилляторов.  
  Например, при значении `1000` `SCINTILLATIONYIELD` будет разделён на 1000.  
  По умолчанию: `1`.

- `-ct, --Trigger-threshold`  
  Нижний порог регистрации энерговыделений в кристалле.  
  По умолчанию: `0` МэВ (значение задаётся в МэВ).

- `-vt, --veto-threshold`  
  Нижний порог регистрации энерговыделений в системе антисовпадений.  
  По умолчанию: `0` МэВ (значение задаётся в МэВ).

- `-oct, --Trigger-optic-threshold`  
  Нижний порог регистрации **оптических фотонов** в кристалле.  
  По умолчанию: `0`.

- `-ovt, --veto-optic-threshold`  
  Нижний порог регистрации оптических фотонов в верхнебоковой системе антисовпадений.  
  По умолчанию: `0`.

- `-obvt, --bottom-veto-threshold`  
  Нижний порог регистрации оптических фотонов в нижней системе антисовпадений.  
  По умолчанию: `0`.

### Геометрия и визуализация

- `-vd, --view-deg`  
  Угол (в градусах), на который разрезается геометрия для просмотра внутренней структуры.  
  По умолчанию: `360`.

- `-vtr, --veto-top-rounded`  
  Радиус кривизны (в мм) для углов верхней части системы антисовпадений.  
  По умолчанию: `0`.

- `-vch, --veto-chamfer-height`  
  Высота (в мм) фаски верхней части системы антисовпадений.  
  По умолчанию: `0`.

- `-vd, --view-deg`  
  Угол (в градусах), на который разрезается геометрия для просмотра внутренней структуры.  
  По умолчанию: `0`.

- `-noUI`  
  Запуск без визуализации.

### Детектор и источник частиц

- `-d, --detector`  
  Тип внутреннего кристалла.  
  По умолчанию: `CsI`.  
  Доступные варианты: `CsI`, `NaI`.

- `-sipm, --Trigger-sipm-config, -csc`  
  Конфигурация SiPM для кристалла.  
  По умолчанию: `12-cross`.  
  Доступные варианты: `2x2`, `8-circle`, `12-cross`, `12-circle`, `12-rhombus`, `13-circle`, `16-cross`.

- `--polished`  
  Делает отражающие поверхности полностью зеркальными.  
  
- `-f, --flux-type`  
  Тип потока частиц.  
  По умолчанию: `Uniform`.  
  Доступные варианты: `Uniform`, `Galactic`, `PLAW`, `COMP`, `Table`, `SEP`.

- `-fd, --f-dir, --flux-dir`  
  Геометрия потока.  
  По умолчанию: `isotropic`.  
  Доступные варианты:  
  `isotropic`, `isotropic_up`, `isotropic_down`,  
  `vertical_up`, `vertical_down`, `horizontal`.

### Дополнительные опции

- `--use-optics`  
  Включает оптическую физику.

- `--save-secondaries`  
  Сохраняет вторичные частицы в отдельный файл.

- `--save-optics`  
  Сохраняет энергию и координаты зарегистрированных фотонов.


### Доступные конфигурации

<div style="display: flex; flex-wrap: wrap; justify-content: center; gap: 20px;">
  <div style="text-align: center;">
    <img src="SiPM_Config_img/Trigger%20SiPM%202x2.png" width="351">
    <br><em>2x2</em>
  </div>
  <div style="text-align: center;">
    <img src="SiPM_Config_img/Trigger%20SiPM%208-circle.png" width="351">
    <br><em>8-circle</em>
  </div>
  <div style="text-align: center;">
    <img src="SiPM_Config_img/Trigger%20SiPM%2012-circle.png" width="351">
    <br><em>12-circle</em>
  </div>
  <div style="text-align: center;">
    <img src="SiPM_Config_img/Trigger%20SiPM%2012-cross.png" width="351">
    <br><em>12-cross</em>
  </div>
  <div style="text-align: center;">
    <img src="SiPM_Config_img/Trigger%20SiPM%2012-rhombus.png" width="351">
    <br><em>12-rhombus</em>
  </div>
  <div style="text-align: center;">
    <img src="SiPM_Config_img/Trigger%20SiPM%2013-circle.png" width="351">
    <br><em>13-circle</em>
  </div>
  <div style="text-align: center;">
    <img src="SiPM_Config_img/Trigger%20SiPM%2016-cross.png" width="351">
    <br><em>16-cross</em>
  </div>
</div>