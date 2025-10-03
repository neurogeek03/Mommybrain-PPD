echo "subdir,missing" > missing_report.csv

for d in */; do
    present=$(find "$d" -maxdepth 1 -type f -name '*.csv' \
              | sed -E 's/.*_subset_([0-9]+)_RCTD_results\.csv/\1/' \
              | sort -n | uniq)
    # Ensure both are sorted the same way
    missing=$(comm -23 <(seq 1 150 | sort -n) <(echo "$present" | sort -n))
    if [ -z "$missing" ]; then
        echo "${d%/},None" >> missing_report.csv
    else
        # Turn spaces into commas for cleaner CSV
        echo "${d%/},\"$(echo $missing | tr ' ' ',')\"" >> missing_report.csv
    fi
done