---
name: biaxnm-inputs
description: BiaxNM 입력 TOML 의 [code], [analysis], [loads] 조각을 만든다. 설계기준, 재료 강도, 횡철근 형식, 하중(Pu, Mux, Muy)을 자연어로 받아 정리한다. 단면 조각(biaxnm-section)보다 먼저 /biaxnm skill 이 호출한다.
---

# 설계기준, 재료, 하중 조각 만들기

결과물은 `[code]`, `[analysis]`, `[loads.이름]` TOML 조각이다. 단면 형상과 배근은 `biaxnm-section` 의 일이다. 단위는 응력 MPa, **하중만 kN, kN·m** 이다.

묻는 것은 AskUserQuestion 으로 하고, 입력에 없는 값은 추측하지 않는다. 호출한 쪽이 넘긴 **이미 아는 값**은 다시 묻지 않고, **부족한 항목**만 묻는다. 아는 값은 표로 보여 주기만 하고 응답을 기다리지 않는다. 부족한 항목이 없으면 질문 없이 바로 조각을 낸다.

이 skill 은 단면보다 먼저 실행되므로 형강·PS강봉이 있는지 아직 모른다. `Fy`(형강 항복강도)와 `fpy`(긴장재 항복강도)는 여기서 묻지 않는다. `biaxnm-section` 이 형강이나 PS강봉을 확인하고 받아서 /biaxnm 이 `[code]` 에 더한다.

## 0. 스키마 읽기

이 문서 맨 아래 "부록. 입력 파일 스키마"의 `[code]`, `[analysis]`, `[loads]` 부분을 읽는다. 설계기준별 키(ACI 는 `fc`, KDS·EC2 는 `fck` 등)는 거기서 확인한다. 키 이름과 기본값을 기억에 의존하지 않는다.

## 1. 설계기준

`type` 을 묻는다.

| type | 기준 | 비고 |
|---|---|---|
| `aci318-25` | ACI 318-25 | `fc`, `fy`, `transverse` |
| `kds14-22` | KDS 14 20 20:2022 강도설계법 | `fck`(90 이하), `fy`, `transverse`. KDS 24 14 20 교량 극한강도설계법도 같음 |
| `kds14-22-app` | KDS 14 20 20:2022 부록(재료계수) | `fck`, `fy`, `phi_c`, `phi_s`, `stress_block` |
| `kds24-25` | KDS 24 14 21 교량 한계상태설계법 | `fck`, `fy`(600 이하), `combination` |
| `ec2-04-1-1` | EN 1992-1-1:2004 (Eurocode 2 건축, 1세대) + National Annex | `fck`(12~90), `fy`, `na`(`EN` 기본, `UK`, `DE`, `FR`, `ES`, `custom`) |
| `ec2-05-2` | EN 1992-2:2005 (Eurocode 2 교량, 1세대) + National Annex | `fck`(EN 30~70), `fy`, `na`(`EN` 기본, `UK`, `FR`, `ES`, `CZ`, `custom`). αcc 권장값 0.85, 철근 등급 B·C |

교량 한계상태설계법이면 `kds24-25`, 건축 콘크리트구조 강도설계법이면 `kds14-22`, Eurocode 2 건축(EN 1992-1-1)이면 `ec2-04-1-1`, Eurocode 2 교량(EN 1992-2)이면 `ec2-05-2` 가 맞는지 사용자 표현으로 확인한다. Eurocode 2 라고만 하고 건축인지 교량인지 모호하면 묻는다. 모호하면 선택지로 묻는다. 2세대 Eurocode(EN 1992-1-1:2023)는 지원하지 않는다고 알린다.

## 2. 재료와 횡철근

- 필수: 콘크리트 강도(`fc` 또는 `fck`), 철근 항복강도 `fy`
- **프로그램이 받는 범위를 먼저 확인한다.** 값을 받으면 TOML 을 만들기 전에 아래 범위를 검사하고, 벗어나면 그 값을 쓰지 말고 사용자에게 알린 뒤 어떻게 할지(값 변경, 다른 설계기준 등) 묻는다. 사용자가 새 값을 주기 전에는 임의로 바꾸지 않는다.
  - `kds14-22`, `kds14-22-app`: `fck` 90 MPa 이하 (초과하면 프로그램이 `fck above 90 MPa ... not supported` 오류를 낸다)
  - `kds24-25`: `fy` 600 MPa 이하 (초과하면 `fy above 600 MPa is not valid` 오류를 낸다)
  - `ec2-04-1-1`: `fck` 12~90 MPa, `fy` 600 MPa 이하 (`na` 가 `DE`, `FR`, `ES` 이면 500 MPa 이하). 벗어나면 `ValueError` 를 낸다
  - `ec2-05-2`: `fck` 는 `na` 별로 `EN` 30~70, `UK` 25~70, `FR`·`CZ` 20~90, `ES` 25~90 MPa. `fy` 는 600 MPa 이하(`FR`, `ES` 는 500). 철근 등급 `ductility` 는 B 또는 C
- `Fy`, `fpy`, `Ep` 는 이 단계에서 넣지 않는다(위 참고).
- `aci318-25`, `kds14-22` 는 `transverse` 를 묻는다(`tied` 또는 `spiral`). 강도감소계수가 달라진다. 단면 조각의 횡철근 형식(띠철근, 나선철근)과 같아야 한다.
- `kds14-22-app`: `stress_block`(`rectangular` 기본 또는 `parabola`)은 사용자가 말했을 때만 넣는다. `phi_c`, `phi_s` 도 마찬가지이다.
- `kds24-25`: `combination`(`ultimate` 기본, `extreme`, `service`, `fatigue`)은 사용자가 말했을 때만 넣는다.
- `ec2-04-1-1`, `ec2-05-2`: `na`(National Annex)를 묻는다. 사용자가 국가를 말하지 않으면 묻지 말고 `EN`(권장값) 을 쓰고 ASSUMPTIONS 에 적는다. 프리셋에 없는 국가는 `na = "custom"` 에 `gamma_c`, `gamma_s`, `alpha_cc` 를 사용자에게 받아 넣는다. `situation`, `stress_block`, `steel_branch`, `ductility`, `eps_ud`, `gamma_c`, `gamma_s`, `alpha_cc` 는 사용자가 말했을 때만 넣는다. 강도감소계수 phi 는 없다(재료계수로 해석).
- 요청에 없는 값(`Es`, `Ep`, `eps_cu`, `gamma_c` 등)은 쓰지 않고 프로그램 기본값에 맡긴다.
- 단위가 MPa 가 아니면(kgf/cm², ksi 등) 환산하고 환산식을 ASSUMPTIONS 에 적는다.

## 3. 하중

- 하중이 있으면 이름(`L1`, `L2`, ...), `Pu`(kN, 압축이 +), `Mux`(x축 둘레, kN·m), `Muy`(y축 둘레, kN·m) 를 받는다. 여러 개 가능하다.
- 단위가 tf, tf·m 이면 환산한다(1 tf = 9.80665 kN). 환산식을 ASSUMPTIONS 에 적는다.
- 인장이면 `Pu` 가 음수이다. 부호를 사용자가 어떻게 썼는지 확인한다.
- 모멘트 기준점은 기본으로 콘크리트 총단면(구멍 제외)의 도심이다. 입력하는 `Mux`, `Muy` 가 같은 점에 대한 값이어야 한다는 것을 알리고, 다른 점이면 `reference` 를 묻는다.
- 하중 없이 곡선만 원하면 `lambda_u_deg`(모멘트 방향 `atan2(Muy, Mux)`, 도) 를 묻는다.

## 4. 조각 출력

```toml
[code]
type = "aci318-25"
fc = 30
fy = 420
transverse = "tied"

[analysis]
n_c = 80
reference = "geometric"

[loads.L1]
Pu = 1800
Mux = 250
Muy = 145
```

- `[analysis]` 는 `n_c = 80`, `reference = "geometric"` 만 쓴다. 사용자가 다른 해석 옵션을 말했을 때만 더 넣는다. 곡선만 원하면 여기에 `lambda_u_deg` 를 넣는다.
- `[output]` 은 쓰지 않는다. 기본 출력 위치가 `out/<파일 이름>/` 이다.
- 요청에 있던 값과 직접 환산한 값을 구분한 ASSUMPTIONS 를 함께 낸다.

## 부록. 입력 파일 스키마 ([code], [analysis], [loads], [output])

입력 파일은 TOML(또는 JSON) 하나이고 표는 `[code]`, `[section]`, `[analysis]`, `[loads.이름]`, `[output]` 이다. 이 문서에 없는 키를 쓰면 프로그램이 오류로 알려 준다(오타 방지). 파일 안의 상대 경로(`csv`, `output.dir`)는 그 파일이 있는 폴더 기준이다. `[section]` 의 키는 `biaxnm-section` 문서에 있다.

### 단위와 부호

- 길이 mm, 면적 mm², 응력과 탄성계수 MPa, 변형률 무차원, 각도 deg. **하중만 kN, kN·m** 이다.
- 축력은 압축이 +, 인장이 -. `Mux` 는 x축 둘레, `Muy` 는 y축 둘레 모멘트이고 방향은 `atan2(Muy, Mux)` 이다.
- 좌표는 전역 (x, y) mm 이고 모든 좌표가 같은 원점을 쓴다.

### [code] 공통 키

| 키 | 의미 | 기본 |
|---|---|---|
| `fy` | 철근 항복강도 [MPa] (필수) | |
| `Es` | 철근 탄성계수 | 200000 |
| `Fy`, `Es_steel` | 형강 항복강도 (형강이 있으면 필수), 탄성계수 | -, 200000 |
| `fpy`, `Ep` | 긴장재 항복강도 (긴장재가 있으면 필수), 탄성계수 | -, 195000 |

### type = "aci318-25"

`fc`(필수, 콘크리트 압축강도 [MPa]), `fy`(필수), 위 공통 키, 그리고:

| 키 | 의미 | 기본 |
|---|---|---|
| `transverse` | `"tied"`(phi_c 0.65) 또는 `"spiral"`(phi_c 0.75) | `tied` |
| `phi_strain` | phi 를 정하는 인장변형률 εt 의 기준. `"all"` 은 철근, 형강, 긴장재 중 가장 많이 늘어난 것, `"rebar"` 는 철근만(없으면 `all`) | `all` |
| `eps_cu` | 극한 압축변형률 | 0.003 |
| `gamma_c`, `gamma_s` | 콘크리트, 철근·형강·긴장재 재료계수(강도를 이 값으로 나눔) | 1.0 |

ACI 318-25 는 phi 에 상한을 둔다(21.2.2.3): 0.1·fc·Ag ≤ Pn ≤ Pn,bal 구간에서 phi 는 (Pn=0.1·fc·Ag 에서 0.90, Pn=Pn,bal 에서 phi_c)를 잇는 직선을 넘을 수 없다. Pn,bal 은 εt 가 항복변형률에 이르는 점의 Pn 이다. 원문은 비긴장재 부재에만 적용하지만 긴장재가 있는 단면에도 보수적으로 적용한다.

### type = "kds14-22" (KDS 14 20 20:2022 강도설계법)

`fck`(필수, 90 이하), `fy`(필수), 공통 키, 그리고 `transverse`(`tied` phi_c 0.65 / `spiral` 0.70), `stress_block`(`"rectangular"` 기본, 등가 직사각형 블록 4.1.1(8) / `"parabola"` 포물선-직선 4.1.1(7)), `phi_max`(phi 상한. 예: 0.75, 프리텐션 부재 끝부분 단면을 보수적으로 검토, KDS 14 20 10 4.2.3(2)). phi 는 철근의 순인장변형률 εt 로 정하고(긴장재가 있으면 최외단 긴장재) 인장지배 0.85, 압축지배 0.65/0.70 이다. εcu, η, β1 은 fck 에서 정해진다(εcu 0.0033, 40 MPa 초과 시 10 MPa 마다 0.0001 감소). KDS 24 14 20(교량 극한강도설계법)도 같은 계산이다.

### type = "kds14-22-app" 또는 "kds24-25" (재료계수 적용)

`fck`(필수), `fy`(필수), 공통 키, 그리고 `phi_c`(콘크리트 재료계수, 기본 0.65), `phi_s`(철근·형강·긴장재 재료계수, 기본 0.90). 콘크리트는 최대 응력 0.85·phi_c·fck 로 해석하고 결과를 그대로 설계강도로 쓴다(phi 를 다시 곱하지 않고, 최대 축강도 0.80/0.85 상한과 최소 편심은 적용하지 않는다).

- `kds14-22-app`: `stress_block`(`rectangular` 기본 / `parabola`). 철근 응력 = phi_s·min(Es·ε, fy). 설계축강도 Pdo 는 부록 3.2(1)의 식 (3-2), (3-3)으로 따로 계산해 상한으로 쓴다(Pnd ≤ Pdo). 긴장재가 있는 단면은 식이 없으므로 해석의 순수압축 값을 쓴다.
- `kds24-25`: `combination` = `"ultimate"`(기본, 0.65/0.90) | `"extreme"` | `"service"` | `"fatigue"`(1.0/1.0). `phi_c`, `phi_s` 를 직접 주면 `combination` 값을 덮어쓴다(품질관리 충분 시 증가 계수 등). 콘크리트는 포물선-직선 곡선만, 철근 응력 = min(Es·ε, phi_s·fy), `fy` 는 600 MPa 이하.

### type = "ec2-04-1-1" (EN 1992-1-1:2004) 또는 "ec2-05-2" (EN 1992-2:2005) + National Annex

두 type 은 해석이 같고 NA 값이 다르다. `fck`(필수, `ec2-04-1-1` 은 12~90, `ec2-05-2` 는 NA 별 범위), `fy`(필수, NA 상한 이하), 공통 키(`fpy` 는 fp0.1k 로 읽고 gamma_s 로 나눈다. `Fy` 는 gamma_a 로 나눈다), 그리고:

| 키 | 의미 | 기본 |
|---|---|---|
| `na` | `ec2-04-1-1`: `"EN"`(권장값), `"UK"`, `"DE"`, `"FR"`, `"ES"`. `ec2-05-2`: `"EN"`, `"UK"`, `"FR"`, `"ES"`, `"CZ"`. 공통으로 `"custom"` (`gamma_c`, `gamma_s`, `alpha_cc` 를 모두 줘야 함) | `EN` |
| `situation` | `"persistent"`(지속/일시) 또는 `"accidental"` | `persistent` |
| `gamma_c`, `gamma_s`, `alpha_cc` | NA 값을 덮어씀. `alpha_cc` 는 0.8~1.0 | NA 값 |
| `stress_block` | `"parabola"`(포물선-직선 3.1.7(1)) 또는 `"rectangular"`(등가 직사각형 3.1.7(3)) | `parabola` |
| `steel_branch` | `"horizontal"`(변형률 제한 없음) 또는 `"inclined"`(경사 가지, 변형률 한계 `eps_ud`, pivot A) | `horizontal` |
| `ductility` | `"A"`, `"B"`, `"C"` (`inclined` 의 k = ft/fy 와 eps_uk). `ec2-05-2` 는 B, C 만 허용 | `B` |
| `eps_ud` | `inclined` 의 변형률 한계 | NA 값(DE 0.025) 또는 0.9·eps_uk |
| `gamma_a` | 형강 부분계수 | 1.0 |

NA 프리셋(Scia 의 NA 요약 자료 기준이라 NA 원문과 대조해야 한다). `ec2-04-1-1`: EN 은 γc 1.5, γs 1.15, αcc 1.0(우발 1.2/1.0). UK 는 αcc 0.85. DE 는 αcc 0.85, 우발 γc 1.3, fy ≤ 500, eps_ud 0.025. FR 은 fy ≤ 500. ES 는 우발 γc 1.3, fy ≤ 500. `ec2-05-2`(EN 1992-2 3.1.6(101)P: αcc 권장값 0.85): EN αcc 0.85, UK 0.85(휨·축력), FR 1.0, ES 1.0, CZ 0.90. γc, γs, fy 상한은 EN 1992-1-1 의 같은 나라 값을 쓴다. fcd = αcc·fck/γc, fyd = fy/γs 로 해석하고 결과를 그대로 설계강도로 쓴다(phi 없음). 변형률 평면은 pivot B(εcu2), pivot C(점 C, εc2), 경사 가지일 때 pivot A(eps_ud)이다. 최소 편심 e0, 세장 효과, 구속 효과는 반영하지 않는다.

### type = "custom" (연구용: 곡선, 기준점, phi 를 직접 지정)

```toml
[code]
type = "custom"
concrete = { csv = "concrete_test.csv", beyond = ["hold", "zero"] }   # 필수
rebar    = { type = "elastic_plastic", E = 200000, fy = 420 }         # 철근이 있으면 필수
steel    = { ... }                                                    # 형강이 있으면 필수
tendon   = { ... }                                                    # 긴장재가 있으면 필수
eps_cu = 0.003        # 극한 압축변형률. 압축 연단의 변형률로 고정된다 (기본 0.003)
phi = 1.0             # 강도감소계수, 상수 (기본 1.0). 축력 상한은 적용하지 않는다
```

phi 를 함수로 주거나 사용자 Pivot, 축력 상한이 필요하면 Python 으로 `CustomCode` 를 쓴다.

### 재료 곡선 (custom 의 concrete/rebar/steel/tendon, 그리고 철근·형강·긴장재의 `law` 에 공통)

압축이 +. 세 가지 형식 중 하나를 쓴다.

1. 점 표: `{ strain = [-0.08, -0.0023, 0.0023, 0.08], stress = [-480, -480, 480, 480] }`
2. CSV: `{ csv = "파일.csv", delimiter = "," }` (열 2개: 변형률, 응력. 머리글과 `#` 주석 행은 건너뜀)
3. 탄소성: `{ type = "elastic_plastic", E = 200000, fy = 420 }`

점 표와 CSV 에는 표 바깥 변형률의 거동 `beyond` 를 줄 수 있다: `"hold"` 는 양쪽 모두 끝 응력 유지, `"zero"` 는 양쪽 모두 응력 0(파단, 압괴), `["hold", "zero"]` 는 [인장 쪽(작은 변형률), 압축 쪽(큰 변형률)] 순서이다. 점 표에는 `yield_strain = 0.0021` 을 줄 수 있다(phi 계산에 쓰는 항복변형률. 없으면 설계기준 값 사용).

### [analysis] 해석 옵션 (전부 선택)

`[loads]` 가 있으면 곡선은 각 하중의 모멘트 방향으로 자동 계산되므로 `lambda_u_deg` 는 필요 없다. `[loads]` 가 없으면 `lambda_u_deg` 를 줘야 하고 그 방향의 곡선만 계산한다.

| 키 | 의미 | 기본 |
|---|---|---|
| `lambda_u_deg` | 모멘트 방향 atan2(Muy, Mux) [deg] | - |
| `n_c` | 중립축 깊이 c 의 개수. 2 이상의 정수 | 80 |
| `c_min_ratio` | 가장 얕은 c (단면 외접 대각선 D 의 배수) | 0.02 |
| `c_max_ratio` | 가장 깊은 c (D 의 배수, `c_min_ratio` 보다 커야 함) | 10 |
| `reference` | 모멘트 기준점. 입력하는 Mux, Muy 도 같은 점에 대한 값이어야 한다. `"geometric"`: 콘크리트 총단면(구멍 제외)의 도심, 철근·형강 제외. `"plastic"`: 소성 중심(순압축 때 철근, 형강, 긴장재까지 포함한 합력점). `[x, y]`: 임의 좌표 [mm] | `geometric` |
| `tol_deg` | λn = λu 를 맞추는 각도 허용오차 [deg] | 1e-6 |
| `scan_deg` | 해를 못 찾을 때 전 구간 탐색 간격 [deg] | 5 |
| `full_scan` | true 면 매번 0~360도를 전부 탐색(느림, 검증용) | false |

### [loads.이름] 하중 (선택, 여러 개 가능, 단위 kN, kN·m)

표 이름이 출력 파일 이름에 쓰인다(`<basename>_<이름>.png/.csv`). 하중마다 DCR(요구/용량 비)을 구하고 DCR ≤ 1 이면 OK 이다. 키는 `Pu`(소요 축력, 압축 +), `Mux`(x축 둘레), `Muy`(y축 둘레).

### [output] 출력 (선택)

결과 폴더 기본값은 `out/<입력 파일 이름>/` (실행한 폴더 기준, 없으면 만듦). 명령줄 `--out 폴더` 를 주면 그 폴더에 저장하고, `--no-save` 는 출력만 한다. `dir`(기본 폴더 대신 쓸 위치, 상대 경로는 입력 파일 위치 기준)와 `basename`(파일 이름 앞부분, 기본은 입력 파일 이름)을 줄 수 있다.
