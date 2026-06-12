#!/usr/bin/env bash
set -e
cd "$(dirname "$0")/scripts"

echo "=== 01. Legacy cluster genetic analysis (1_clustering) ==="
python3.11 01_legacy_clusters.py

echo ""
echo "=== 02. GMM cluster genetic analysis ==="
python3.11 02_gmm_clusters.py

echo ""
echo "=== 03. Legacy cluster analysis (NT restricted to C1) ==="
python3.11 03_legacy_clusters_nt_c1.py

echo ""
echo "=== 04. GMM cluster analysis (NT restricted to C1) ==="
python3.11 04_gmm_clusters_nt_c1.py

echo ""
echo "=== 05. hg38 Legacy cluster genetic analysis ==="
python3.11 05_hg38_legacy_clusters.py

echo ""
echo "=== 06. hg38 GMM cluster genetic analysis ==="
python3.11 06_hg38_gmm_clusters.py

echo ""
echo "=== 07. hg38 Legacy cluster analysis (NT restricted to C1) ==="
python3.11 07_hg38_legacy_clusters_nt_c1.py

echo ""
echo "=== 08. hg38 GMM cluster analysis (NT restricted to C1) ==="
python3.11 08_hg38_gmm_clusters_nt_c1.py

echo ""
echo "=== 09. SPARK LoF carrier analysis (GMM clusters) ==="
python3.11 09_spark_lof_clusters.py

echo ""
echo "=== 10. SPARK LoF carrier analysis (Manual clusters) ==="
python3.11 10_spark_lof_manual_clusters.py

echo ""
echo "=== 11. Manual cluster genetic analysis (hg19) ==="
python3.11 11_manual_clusters.py

echo ""
echo "=== 12. hg38 Manual cluster genetic analysis ==="
python3.11 12_hg38_manual_clusters.py

echo ""
echo "=== 13. SPARK LoF with LEAP-InovAND NT (GMM clusters) ==="
python3.11 13_spark_lof_leapnt_clusters.py

echo ""
echo "=== 14. SPARK LoF with LEAP-InovAND NT (Manual clusters) ==="
python3.11 14_spark_lof_leapnt_manual_clusters.py

echo ""
echo "=== All done ==="
