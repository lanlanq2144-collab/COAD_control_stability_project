# Notes

This file contains study notes related to Python, GitHub, Kaggle, statistics, and biomedical data analysis for this project.

## GitHub Notes

* A repository is a project folder that stores files, code, documentation, and version history.
* A commit saves a version of the project at a specific point in time.
* `README.md` is the main introduction page for the project.
* `work_log.md` is used to document time spent and project progress.
* `project_plan.md` is used to organize the research question, data sources, methods, and timeline.

## Python Notes

* Python will be used for data cleaning, analysis, visualization, and machine learning.
* `pandas` is commonly used to work with table-like data.
* `NumPy` is commonly used for numerical calculations.
* `matplotlib` is commonly used to create graphs and figures.
* `scikit-learn` is commonly used for machine learning models and feature selection.

## Research Notes

* Colon adenocarcinoma is commonly abbreviated as COAD.
* Tumor-adjacent normal tissue may look normal under a microscope, but it may not be molecularly identical to healthy normal tissue.
* This project focuses on whether biomarker candidates selected by machine learning remain stable when the normal control group changes.
* The planned comparison is between tumor-adjacent normal tissue and healthy GTEx colon tissue.

## Data Source Notes

* TCGA can provide tumor and tumor-adjacent normal samples.
* GTEx can provide healthy normal tissue samples.
* UCSC Xena can be used to access and download TCGA and GTEx gene expression data.

## Questions to Explore

* How different are tumor-adjacent normal tissue and healthy normal tissue at the gene expression level?
* Do machine learning-selected biomarker candidates change when the normal control group changes?
* Which genes remain stable across different control group comparisons?
* How should batch effects between TCGA and GTEx be handled?

Progress on Jun 8(but in Korean)

-연구 

-1단계: QC(Quality Control) 품질 확인 (데이터가 신뢰 가능한지 품질을 확인하는 단계)

1-1: Expression matrix 생성
: UCSC Xena에서 엑셀 자료를 다운받고 거기에서 대장암 관련 자료만 엑셀로 추출했다.(19131개 중 637개만 * 60498가지의 유전자), chat gpt에게 파일을 만들어달라고 요청했고 gpt에서 코드를 줬다.

1-2: 전체 PCA
: 배치효과를 확인하기위해 전체 PCA를 찍었고 60498개의 유전자를 2D로 정리하는거라고는 하는데 아직 Linear Algebra를 배우지 않은상태라 그런가 아직 감은 잘 잡히지 않기에 여름학기에 들으면 이해가 될 것같다고 생각한다.
PC1: 24.9%라는건 이 축 하나가 전체 정보의 24.9프로를 설명한다는것이다. 
4개 그룹을 한번에 찍었어요. 결과에서 GTEx가 왼쪽에 따로 뭉치고 TCGA가 오른쪽에 있었어요. 여기서 배치효과가 있다는 걸 처음 확인했어요.

네 맞아요. 점이 멀리 떨어져 있을수록 그 그룹들 사이의 차이가 크다는 뜻이에요. 근데 그 차이가 생물학적 차이인지 배치효과인지는 PCA만으로는 구분이 안 돼요.

배치효과가 뭐냐
같은 조직인데 측정 환경이 달라서 생기는 인위적인 차이예요.
예를 들어 같은 사람 혈액을 A병원 장비로 재면 혈당 100, B병원 장비로 재면 혈당 115가 나온다고 해보자. 혈당이 실제로 다른 게 아니라 장비가 달라서 숫자가 다르게 나오는 거예요. 이게 배치효과예요.

우리 연구에서
TCGA  →  A 기관, A 장비, A 방식으로 측정
GTEx  →  B 기관, B 장비, B 방식으로 측정
둘 다 정상 대장 조직인데 측정 환경이 달라서 PCA에서 멀리 떨어진 거예요. 이게 배치효과예요.

그래서 NAT vs GTEx PCA에서
NAT랑 GTEx가 멀리 떨어져 있었는데, 둘 다 정상 조직이니까 이 분리가 생물학적 차이보다 배치효과일 가능성이 높다고 판단한 거예요.
이해됐어요?

(PCA(Principal Component Analysis, 주성분 분석)는 데이터의 정보 손실을 최소화하면서 고차원 데이터를 저차원 공간으로 축소하는 대표적인 차원 축소(Dimensionality Reduction) 기법)라고 인터넷에서 말하고있다.

1-3: TCGA내부 배치효과 확인(NAT vs Tumor)
: NAT, Tumor둘 다 분리가 잘 되어있음=암조직과 정상조직은 당연히 분리되어있어야 함=배치효과 없음

1-4: NAT, GTEx 분리
: NAT, GTEx 둘 다 정상조직인데 분리됨= 배치효과 의심

1-5: Combat을 통해 배치효과를 최대한 제거했음(NAT, GTEx에서)
: 그러나 몇 개 점들은 여전히 멀리 치우쳐있으나, 이유는 배치효과 제거가 덜 되거나, 암환자의 몸속이라 정상인들과 무언가 다르다고 추측했음.
또한 Tumor는 왜 combat을 처리하지 않았냐면, 암과 정상의 생물학적 차이까지 지워버리면 힘들어지기때문.


-2단계: DEG(Differentially Expressed Genes, 차등 발현 유전자), 암조직, 정상조직에서 발현이 다른 유전자 찾기

2-1: 데이터 분리(Tumor vs NAT, Tumor vs GTEx)
: Tumor vs NAT= 둘 다 TCGA데이터이기때문에 Batch효과가 없어서 원본 데이터로
: Tumor vs GTEx= 암은 TCGA, 정상은 GTEx에서 가져왔고 Batch효과를 확인했음=Combat보정된 데이터를 써야 함
>claude가 이 부분이 보고서에서 공격받기 쉬운 방법론이라고 했는데 이유는
: Combat을 NAT + GTEx에만 적용했기 때문이라고.
: 배치효과를 제거해서 GTEx가 TCGA스케일에 맞춰졌고, Tumor는 ComBat에 포함 안 되었고
??: 보정된 GTEx랑 안 된 Tumor를 비교하는게 맞나? 라고 물으면
>왜 Tumor를 제거했는가: Tumor를 ComBat에 넣으면, TCGA, GTEx의 배치효과를 제거해서 Tumor와 정상의 진짜 생물학적 차이까지 지워버릴 수 있기때문에 어쩔 수 없이 뺐다.
== 완벽한 방법이 없다. 방법론적 질문을 받을 수 있기에 limitation에 솔직히 쓰는게 좋음.

2-2: 통계 검정
: Tumor샘플이 288개이므로 288명의 통계를 검정했다.
log2FC가 뭔지
FC는 Fold Change, 몇 배 차이나는지예요.
예를 들어:
유전자 A
암 환자 평균   → 100
정상 평균      → 10
→ 10배 차이
근데 10배를 그냥 쓰면 숫자가 너무 커지고 다루기 불편해요. 그래서 log2를 씌워요.

log2(10배) = 약 3.3
log2(2배)  = 1
log2(0.5배) = -1  ← 정상에서 높은 것
그래서 log2FC가 양수면 암에서 높고, 음수면 정상에서 높은 거예요.

p-value가 뭔지
이 차이가 우연히 생겼을 확률이에요.
p-value = 0.001
→ 이 차이가 우연일 확률이 0.1%
p-value = 0.9
→ 이 차이가 우연일 확률이 90%

2-3: 필터링
: 기준을 세웠고, 기준을 통과한 유전자만 DEG로 선택했다. (60,498개는 너무 많고 신뢰가 어렵기때문에 두 가지 동시 통과한 유전자만)
ㅣlog2FCㅣ> 1, FDR < 0.05(p-value를 쓰면 p-value가 낮게 나오는게 수천개씩 나옴. 그래서 FDR(False Discovery Rate)은 수만 개를 동시에 볼 때 보정한 값이므로 p-value보다 더 엄격함, 실제로는 차이 없는 게 몇 %냐를 조절함

2-4: Volcano Plot(DEG결과 그래프 시각화)
:PCA를 할때는 GPT로 했지만 이때부턴 Claude Code를 결제하여 그래프를 만들어달라고 요청함
log2FC = log2(암 평균 / 정상 평균)이기에 x축이 오른쪽으로가면 암에서 더 많이 켜져있고, x축에서 왼쪽으로가면 정상에서 더 많이 켜져있다는것이 된다. 
y축은 이 차이가 우연이 아닐 가능성이 높다는것을 나타내고. FDR이 작을수록 정확하기에 -스케일을 적용한다.

3단계: ML(Machine Learning)
: DEG에서 나온 유전자 수천 개 중에서 “암 판단에 진짜 핵심인 유전자”를 추리기 위함.

3-1: 데이터 준비(모두 ML에 넣으면 문제가 발생하기에 데이터를 추리는것이 중요함)
: 그래서 상위 500개를 골랐음(FDR이 가장 낮은)
라벨을 붙였는데 
Tumor Sample = 1
Normal Sample =0
라고 붙이면 ML 모델이 “이 샘플이 암이다 아니다”를 알려주는 것이고, 모델이 이를 학습한다고함.
표를 가져와 1-1의 Expressioin Matrix에서 500개 행만 추출하여 라벨열 하나를 추가함, Label= 1 or 0
(Claude code가 Pandas로 처리)

3-2: LASSO(최소한의 유전자로 암과 정상을 구분하는 방법)
500개의 유전자중 중요하지 않은 유전자의 계수를 0으로 만듦, 이로인해 500개중에서 1~3개만 남아버림, 기준을 덜 엄격하게했을땐 19~40개가 나옴.
1~3개짜리로 하면 “이 유전자들이 있으면 암인지 정상인지 구분 가능하다”는 증명 but 바이오마커 후보로는 너무 적기에 19~40개 중 Random Forest결과와 겹치는 것들만 최종 후보로 선택.

어떻게 작동하냐
500개 유전자가 있는데 전부 다 중요한 건 아니에요. LASSO는 중요하지 않은 유전자의 계수를 0으로 만들어버려요.


유전자1  →  계수 2.3  (중요함, 살아남음)
유전자2  →  계수 0    (중요하지 않음, 제거됨)
유전자3  →  계수 -1.5 (중요함, 살아남음)
유전자4  →  계수 0    (제거됨)



3-3: Random Forest(여러 개의 결정 트리를 만들어서 투표하는 방법)
결정 트리가 뭔지 먼저
하나의 규칙 묶음이에요.


예시:
유전자A > 3이면 → 암
유전자A < 3이면 → 유전자B 확인
  유전자B > 5이면 → 암
  유전자B < 5이면 → 정상
이런 if/else 규칙이 트리 모양으로 생겼어요.

Random Forest는
이런 트리를 300개 만들어요. 근데 각 트리마다 쓰는 유전자랑 샘플을 랜덤하게 조금씩 다르게 해요.


트리1: "유전자A > 3이면 암"
트리2: "유전자B > 7이면 암"
트리3: "유전자A > 2 이고 유전자C > 4면 암"
...
300개 트리가 각자 판단
→ 다수결로 최종 결정

중요도가 어떻게 계산되냐
300개 트리에서 어떤 유전자가 자주 쓰였는지 점수 매겨요.


유전자A → 280개 트리에서 쓰임 → 중요도 높음
유전자B → 10개 트리에서 쓰임  → 중요도 낮음

LASSO랑 차이


LASSO        →  유전자들을 더하거나 빼서 판단 (선형)
Random Forest →  유전자 조합 패턴으로 판단 (비선형)
두 방법이 다른 시각으로 보는 거라서 둘 다 동의한 유전자가 더 신뢰할 수 있어요.

그렇게 LASSO 느슨한 버전 19~40개, Random Forest 상위 50개를 합쳐서 교집합을 만들었는데, 둘 다 중요하다고 동의한 유전자가 7~14개 나옴(두 방법이 서로 다른 시각으로 봤는데 둘 다 중요하다고 한 유전자이기에 더 신뢰가능함) 

3-4: 모델 성능 평가(ML 모델이 암, 정상을 얼마나 잘 구분하는지 확인하는 단계)
Train set(모델이 학습하는 데이터), Test set(모델이 본 적 없는 데이터로 테스트)
모델이 학습에 쓴 데이터는 당연히 잘 맞추기에 본 적 없는 데이터로 테스트하여 진짜 성능을 확인

AUC란: 0.5면 동전던지기 수준, 1.0이어야 완벽한 구분. 4개 모두 AUC = 1.000이 나옴{엥?}

왜 1.000이 나왔나?: DEG상위 500개를 입력으로 썼기에, 이미 암으랑 정상에서 가장 극단적으로 차이나는 유전자들이기에 모델이 쉽게 구분 가능함.

문제: Feature-Selection Leakage 가능성
전체 데이터에서 DEG 500개 먼저 뽑음
→ 그 500개로 train/test 나눠서 평가
→ test 데이터 정보가 이미 DEG 선택에 
  새어들어간 것
→ 성능이 과장됐을 수 있음
그래서 보고서 limitation에 이걸 써야 해요.(라는데 데이터 초기화하고 다시 하면 되는거 아닌가)
라는 의문이 들어서 claude에게 물어보고 nested cross-validation를 실행하여 정확도를 올리기로 하고, claude code로 실행했으나 결과 동일 
leakage수정후에도 AUC 1.000
1. 방법론적으로 깨끗하다는 증거
2. 생물학적 신호가 진짜로 강해서 어떻게 테스트해도 완벽하게 구분
보고서 limitation에도 원래 AUC leakage가능성 있음에서 그 자리에 leakage 검증완료했으며 결과는 동일함

3-5 바이오마커 후보 추출(LASSO, Random Forest 결과를 합쳐서 최종 후보를 뽑는 단계)
왜 두개를 합치냐?: LASSO만 쓰면 놓친 유전자 있을 수 있음, Random Forest만 쓰면 RF가 놓친 유전자가 있을수있음, 둘 다 동의했으면 두 방법이 서로 다른 시각으로 봤는데 둘 다 중요하다고 한것 > 더 신뢰 가능
교집합 뽑기


LASSO (느슨한 버전) 19~40개
∩
Random Forest 상위 50개
= 7~14개 / 비교군

4개 비교군에서 각각 나온 후보


Tumor vs NAT        →  후보 목록 A
Tumor vs GTEx 전체  →  후보 목록 B
Tumor vs Transverse →  후보 목록 C
Tumor vs Sigmoid    →  후보 목록 D
(7~14개가 후보)

4번의 ML에서 LASSO에서 나온것 and RF 상위 50개의 교집합만 남겨서 최종 후보 목록을 만듦

4단계: 바이오마커 안정성 비교(3단계에서 4개 비교군 각각에서 후보 목록이 나왔고, 4단계에서는 그 4개 목록을 비교)

4-1: 벤 다이어그램
venn4_biomarker_overlap.png(4개 그룹 전체 한번에)
NAT(분홍)      → 독자적으로만 5개
GTEx_All(파랑) → 독자적으로만 2개
Transverse(초록) → 독자적으로만 5개
Sigmoid(보라)  → 독자적으로만 7개
가운데 겹치는 부분 → 거의 없음
(4개 전부에 동시에 들어있는 유전자가 없음, 그러나!! 4개말고 3개나 2개에서는 겹침)
ex) NAT에서 5개라는건 후보 유전자들중 어디에도 없는 유전자가 5개 있다는 것.
= 정상을 NAT로 정의할때만 나오는 것, GTEx로 바꾸면 사라짐

여기에서 venndiagram에서 NAT독자적으로 5개랑, NAT 바이오마커 후보 파일에서는 7개라고 나오는데 왜 숫자가 다른지 갭이 생김
: NAT후보들중 나머지 2개는 다른 비교군들과도 겹침, 4개 모두 겹치는게 없다는것이지만 3개나 2개는 꽤 있음

4점 = 4개 모두에서 나오는것 그러나 없음
3점 = 3개 비교군에서 나오는 것: 그 중에 2개는 가짜 유전자라서 제거(가짜유전자: ENSG00000279473: 데이터베이스 조회했더니 현존하지 않는 삭제된 유전자라고함, NAT대비 log2FC는 +0.37이나 GTEx대비 +11.6이기에 배치효과 artifact일 확률 높음, HMGN1P36: 유사유전자로 진화 과정에서 기능을 잃은 유전자 잔재로 단백질을 만들지 못함, 마찬가지고 NAT대비 차이는 거의 없으나 GTEx대비로만 극단적으로 높음), CDH3, CLDN1, ETV4, KRT80남음

venn_pairwise.png(둘씩 짝지어서) 
NAT ∩ GTEx_All     → 2개
NAT ∩ Transverse   → 1개
NAT ∩ Sigmoid      → 1개
GTEx_All ∩ Transverse → 8개
GTEx_All ∩ Sigmoid    → 6개
Transverse ∩ Sigmoid  → 5개

== NAT랑 GTEx끼리는 1~2개밖에 안 겹침, GTEx끼리는 5~8개 겹침 = 정상을 NAT로 잡느냐 GTEx로 잡느냐에 따라 바이오마커가 완전히 달라진다는 게 이 그림으로 나타남.

5단계: 생물학적 검증(최종 4개 유전자가 진짜 생물학적으로 의미가 있는지 확인하는 단계)

5-1: 유전자 이름 변환(Ensembl ID를 실제 유전자 이름으로 바꿈)
ENSG00000062038  →  CDH3(Cadherin-3, 세포끼리 붙어있게 하는 접착제 역할, 암: 과발현> 세포가 떨어져 나가 다른 곳으로 이동> 전이랑 관련)
ENSG00000175832  →  ETV4(ETS Variant Transcription Factor 4, 유전자 발현을 조절하는 전사인자, 암: 과발현> 다른 유전자들을 잘못 켜거나 끔> 암세포 증식, 전이 촉진)
ENSG00000163347  →  CLDN1(Claudin-1, 세포 사이 틈을 막는 장벽 역할, 암: 과발현되면서 장벽 기능이 망가져서 암세포가 빠져나가기 쉬워짐)
ENSG00000167767  →  KRT80(Keratin 80, 세포골격을 구성하는 단백질, 암: 과발현되면 AKT 신호경로 활성화> 암세포 침습, 전이 촉진)

5-2: 선행연구 확인(PubMed에서 4개 유전자 각각을 대장암이랑 연결한 논문이 있는지 Claude에게 맡김, 또 읽어봄)

CDH3


대장암 바이오마커로 가장 많이 연구됨
→ TCGA 데이터로 대장암 과발현 확인한 논문 다수
→ 혈액에서도 검출 가능 (AUC = 0.90)
→ 전이 환자에서 특히 높음

CLDN1


문헌 지지가 제일 강함
→ 메타분석으로 예후 바이오마커 공식 확인
→ BRAF 돌연변이랑 연관
→ 항암제 내성과도 연결

ETV4


대장암 특이적 역할 확인
→ 정상 대장 조직에서 암으로 진행될 때 핵심 역할
→ 림프절 전이, 생존 단축과 연관

KRT80


대장암에서 직접 연구된 논문 있음
→ AKT 경로 통해 침습, 전이 촉진
→ siRNA로 억제했더니 암세포 증식 감소

제일 중요한 발견
CDH3 + ETV4 + CLDN1 조합이 다른 팀의 ML 연구(Guo et al. 2022)에서도 독립적으로 나왔어요.


다른 팀이 완전히 다른 방법으로 분석했는데
같은 유전자 조합이 나온 거예요
→ 우리 결과가 재현된 것


5-3: Pathway 분석(4개 유전자가 어떤 생물학적 경로랑 관련있는지 확인하는 단계)

Pathway란: 유전자들은 혼자 일하지 않고 여러 유전자가 팀을 이뤄서 특정 기능을 담당하기에 그팀을 Pathway라 함. 
ex) 예시:
면역 반응 pathway → 면역 관련 유전자들 모임
세포 분열 pathway → 세포 증식 관련 유전자들 모임

우리 4개 유전자가 속한 pathway


Epithelial cell differentiation
→ 상피세포 분화 경로
→ 정상 상피세포가 어떻게 만들어지는지

Skin/Epidermis development
→ 상피세포 구조 관련

Tight junction
→ CLDN1이 세포 사이 틈 막는 것

이게 대장암이랑 어떻게 연결되냐


정상 대장 상피세포
→ 분화가 잘 됨
→ 세포끼리 잘 붙어있음
→ 장벽 기능 정상

암이 되면
→ 분화가 망가짐
→ 세포가 떨어져 나감 (전이)
→ 장벽 기능 붕괴
4개 유전자 전부가 이 과정과 관련있어요.
라고함

5-4: 생존 분석(이 유전자들이 높은 환자, 낮은 환자의 생존율이 다른지 확인한 단계)

각 유전자마다: 발현이 높은 환자 그룹, 발현이 낮은 환자 그룹 > 두 그룹의 생존율 비교
CDH3   → p = 0.90  (유의미하지 않음)
CLDN1  → p = 0.61  (유의미하지 않음)
ETV4   → p = 0.63  (유의미하지 않음)
KRT80  → p = 0.21  (유의미하지 않음)
이유가 두 가지예요.


1. 중간값으로 높음/낮음 나눈 게 너무 단순해요
   → 실제로는 Cox regression으로
     다른 변수들(나이, 병기 등)도 
     같이 고려해야 해요

2. 바이오마커가 진단용이지 예후용이 아닐 수 있어요
   → 암인지 아닌지 구분하는 건 잘 하지만
     생존율 예측은 다른 문제


6단계: 보고서

Claude Code가 자동으로 작성한 영어 논문 형식 초안이에요.

구성


Abstract      → 연구 전체 250단어 요약
Introduction  → 연구 배경, 왜 이 연구가 필요한지
Methods       → 1~5단계 분석 방법 전부
Results       → 결과 (그림 참조 포함)
Discussion    → 결과 해석, 의의, 한계
Conclusion    → 최종 정리
References    → 인용 논문 20개
	
	
GPT 검토 결과 8/10
지적받은 것들:


1. 빈 표 있음          → 채워야 함
2. ComBat 방어 약함    → 보강 필요
3. References 형식     → 일부 PMC 번호만 있음
4. AUC leakage        → 이미 해결됨

지금 남은 것


1. 빈 표 채우기
2. ComBat limitation 문장 보강
3. References 형식 통일

Q. 실제 암 발견할때는 저렇게 왜 안함?(소변이라던가 그런거 검사할때 대략 위험군이다 정도는 알려줄수있는거아냐?)


















