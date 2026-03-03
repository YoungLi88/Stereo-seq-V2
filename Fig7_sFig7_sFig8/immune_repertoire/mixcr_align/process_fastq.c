#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

#define MAX_LINE_LENGTH 1000000

// Function to process a single FASTQ record
void process_record(char *header, char *seq, char *plus, char *qual, gzFile output_file) {
    // Remove trailing newline characters
    header[strcspn(header, "\n")] = '\0';
    seq[strcspn(seq, "\n")] = '\0';
    plus[strcspn(plus, "\n")] = '\0';
    qual[strcspn(qual, "\n")] = '\0';

    // Extract UMI and modify header
    char umi[10];
    strncpy(umi, seq, 9);
    umi[9] = '\0';

    char new_header[MAX_LINE_LENGTH];
    snprintf(new_header, MAX_LINE_LENGTH, "%s|||UB:Z:%s", header, umi);

    // Remove the first 9 bases and corresponding quality values
    char new_seq[MAX_LINE_LENGTH];
    char new_qual[MAX_LINE_LENGTH];
    strncpy(new_seq, seq + 9, strlen(seq) - 9);
    new_seq[strlen(seq) - 9] = '\0';
    strncpy(new_qual, qual + 9, strlen(qual) - 9);
    new_qual[strlen(qual) - 9] = '\0';

    // Write the modified record to the output file
    gzprintf(output_file, "%s\n%s\n%s\n%s\n", new_header, new_seq, plus, new_qual);
}

// Function to process the input FASTQ file
void process_fastq_file(const char *input_file_path, const char *output_file_path) {
    // Open the input and output files
    gzFile input_file = gzopen(input_file_path, "r");
    if (input_file == NULL) {
        perror("Error opening input file");
        exit(EXIT_FAILURE);
    }

    gzFile output_file = gzopen(output_file_path, "w");
    if (output_file == NULL) {
        perror("Error opening output file");
        exit(EXIT_FAILURE);
    }

    char header[MAX_LINE_LENGTH];
    char seq[MAX_LINE_LENGTH];
    char plus[MAX_LINE_LENGTH];
    char qual[MAX_LINE_LENGTH];

    // Read and process the input file line by line
    while (1) {
        if (gzgets(input_file, header, MAX_LINE_LENGTH) == NULL) break;
        if (gzgets(input_file, seq, MAX_LINE_LENGTH) == NULL) break;
        if (gzgets(input_file, plus, MAX_LINE_LENGTH) == NULL) break;
        if (gzgets(input_file, qual, MAX_LINE_LENGTH) == NULL) break;

        process_record(header, seq, plus, qual, output_file);
    }

    // Close the input and output files
    gzclose(input_file);
    gzclose(output_file);
}

int main(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <input.fastq.gz> <output.fastq.gz>\n", argv[0]);
        return EXIT_FAILURE;
    }

    const char *input_file_path = argv[1];
    const char *output_file_path = argv[2];

    process_fastq_file(input_file_path, output_file_path);

    return EXIT_SUCCESS;
}
