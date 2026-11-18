# BiaxNM

임의 다각형 단면(RC, SRC, 긴장재 포함)의 이축휨 P-M 상관도와 DCR을 계산한다. 지원하는 설계기준(`[code] type`)은 아래와 같다.

| type | 기준 |
|---|---|
| `aci318-25` | ACI 318-25 |
| `kds14-22` | KDS 14 20 20:2022 강도설계법 |
| `kds14-22-app` | KDS 14 20 20:2022 부록|
| `kds24-25` | KDS 24 14 21 교량 한계상태설계법 |
| `ec2-04-1-1` | EN 1992-1-1:2004 |
| `ec2-05-2` | EN 1992-2:2005  |


> **주의 및 면책**
> - 실무에 쓰기 전에 설계기준 원문과 독립적인 계산으로 반드시 확인해야 한다.
> - 세장 효과, 최소 편심(Eurocode), 철근비 제한, 상세 규정은 검토하지 않는다. 프리텐션 부재 끝부분의 φ 제한(표 21.2.3)과 긴장재 전달길이도 검토하지 않는다.
> - 이 프로그램은 어떠한 보증 없이 "있는 그대로" 제공된다. 입력 오류, 결과의 해석, 설계와 시공을 포함한 **모든 사용과 그 결과에 대한 책임은 전적으로 사용자에게 있으며**, 작성자는 이 프로그램의 사용으로 생긴 어떠한 손해에 대해서도 책임지지 않는다.

## Installation

Python 3.11 이상과 numpy, scipy, matplotlib이 필요하다.

1. 저장소를 내려받고 폴더로 이동한다.

```bash
git clone https://github.com/Memorhee/BiaxNM.git
cd BiaxNM
```

2. conda 환경을 만들고 활성화한다.

```bash
conda create -n biaxnm python=3.11 numpy scipy matplotlib
conda activate biaxnm
```

conda를 쓰지 않으면 venv와 pip로도 된다.

```bash
python -m venv .venv
.venv\Scripts\activate
pip install numpy scipy matplotlib
```

3. 설치를 확인한다.

```bash
python src/main.py examples/rect_aci318-25.toml --no-save
```

## 실행

```bash
python src/main.py examples/rect_aci318-25.toml
```

Windows에서는 `BiaxNM.bat`을 실행하면 `examples/`의 예제 6개를 차례로 계산한다.

### 옵션

| 옵션 | 설명 |
|---|---|
| `--out 폴더` | 결과 저장 위치 지정 |
| `--no-save` | 출력만 하고 파일은 저장하지 않음 |

### 입력 파일

`examples/rect_aci318-25.toml`의 주석에 있다

| 예제 | 내용 |
|---|---|
| `rect_aci318-25` | 사각 RC 기둥 (ACI 318-25, 모든 옵션 설명) |
| `circle_kds14-22` | 속찬 원형, 나선철근 (KDS 14 20 20) |
| `hollow_kds14-22-app` | 중공 원형 (KDS 14 20 20 부록, 재료계수) |
| `track_ec2-05-2` | 속찬 트랙형 (EN 1992-2) |
| `phc_kds24-25` | 중공 PHC 말뚝, 긴장재 (KDS 24 14 21) |
| `src_ec2-04-1-1` | 형강 매입(SRC) (EN 1992-1-1) |


| 구분 | 선택지 |
|---|---|
| 콘크리트 단면 | 사각, 원형, 트랙형, T형, L형, 임의 다각형(구멍 포함) |
| 철골 | H형, 중공 사각형, 임의 다각형 |
| 철근 | 좌표 + `area` / `dia` / 규격 호칭 `size` (KS, ASTM, ASTM_SI, BS), 선 배근, 원 배근 |
| 긴장재 | 면적 + 유효 프리스트레스 `fse` |
| 재료 곡선 | ACI 기본, 점열, CSV, 탄소성 |


### Claude Code Skill (`/biaxnm`)

Claude Code로 이 저장소 폴더를 열고 `/biaxnm` 을 입력하면 된다. 스킬은 `.claude/skills/` 에 들어 있어 따로 설치하지 않는다. 터미널이든 데스크톱 앱이든 같다. Python 환경(위 Installation)은 미리 갖춰야 한다.

#### 예: 중공 사각 기둥

프롬프트 (철근은 철근비만 준다)

```
/biaxnm 외곽 800x1000 중공 사각 기둥, 속 구멍 400x600 (가운데, 벽 두께 200). ACI 318-25, fc 35, fy 420, 띠철근. 철근비 1~1.5% 정도. 하중 Pu 6000kN, Mux 1200, Muy 800
```

AI 질문

| 단계 | 내용 |
|---|---|
| 되물을 것 | 순피복, 띠철근 직경, 철근 규격표, 바깥면과 안쪽면 철근 구성 (중공 단면은 면별로 묻는다), 형강·PS강봉 유무 |
| 답 (가정) | 순피복 40, 띠철근 D10, KS, 바깥면·안쪽면 모두 D19, 없음 |
| 제안 | 중공 단면은 배근 스크립트가 없어서 직접 계산한다. 간격 250 이하 기준으로 바깥면 14개 + 안쪽면 12개 (D19, 총 26개), 철근비 1.33% (목표 1~1.5% 안), 순간격 경고 없음. 같은 계산에서 간격 200 이하는 1.64% 로 범위 밖, 300 이하는 1.13% |

AI의 실제 질문은 다를 수 있다.
