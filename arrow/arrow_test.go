package arrow

import (
	"log"
	"os"
	"sync"
	"testing"

	"github.com/apache/arrow/go/v14/arrow/array"
	"github.com/apache/arrow/go/v14/arrow/ipc"
)

func readAndVeryifyArrowFile(filePath string, fieldNames []string, batchSize int, rows [][]uint16) {
	// Open the file for reading
	file, err := os.Open(filePath)
	if err != nil {
		log.Fatal(err)
	}
	defer file.Close()

	// Create a new IPC reader
	reader, err := ipc.NewFileReader(file)
	if err != nil {
		log.Fatal(err)
	}
	defer reader.Close()
	defer os.Remove(filePath)

	for i := 0; i < reader.NumRecords(); i++ {
		record, err := reader.Record(i)
		if err != nil {
			log.Fatal(err)
		}

		if record.NumCols() != int64(len(fieldNames)) {
			log.Fatal("Number of columns in record does not match number of field names")
		}

		for colIdx, col := range record.Columns() {
			arr, ok := col.(*array.Uint16)
			if !ok {
				log.Fatal("Column is not of type Uint16")
			}

			for rowIdx := 0; rowIdx < batchSize; rowIdx++ {
				originalIdx := batchSize*i + rowIdx

				if arr.Len() == rowIdx {
					if originalIdx == len(rows) {
						break
					}

					log.Fatal("Column length does not match number of rows")
				}

				foundValue := arr.Value(rowIdx)

				expectdValue := rows[originalIdx][colIdx]

				if foundValue != expectdValue {
					log.Fatal("Value in column does not match value in row ", foundValue, expectdValue)
				}
			}
		}
	}
}

func TestArrowWriteRead(t *testing.T) {
	batchSize := 5
	fieldNames := []string{"Sample1", "Sample2", "Sample3"}
	filePath := "test_matrix.feather"
	writer, err := NewArrowWriter(filePath, fieldNames, batchSize)
	if err != nil {
		log.Fatal(err)
	}
	// Example with 100 rows
	rows := make([][]uint16, 100)
	for i := 0; i < 100; i++ {
		rows[i] = []uint16{uint16(i), uint16(i + 1), uint16(i + 2)}
	}

	for _, row := range rows {
		if err := writer.Write(row); err != nil {
			log.Fatal(err)
		}
	}

	if err := writer.Close(); err != nil {
		log.Fatal(err)
	}

	readAndVeryifyArrowFile(filePath, fieldNames, batchSize, rows)
}

func TestArrowWriterConcurrency(t *testing.T) {
	fieldNames := []string{"Field1", "Field2"}
	writer, err := NewArrowWriter("concurrent_output.feather", fieldNames, 10)
	if err != nil {
		t.Fatal(err)
	}
	defer writer.Close()

	var wg sync.WaitGroup
	numGoroutines := 5
	numWritesPerRoutine := 10

	for i := 0; i < numGoroutines; i++ {
		wg.Add(1)
		go func(writer *ArrowWriter, routineID int) {
			defer wg.Done()
			for j := 0; j < numWritesPerRoutine; j++ {
				// rowIndex := "row_" + strconv.Itoa(routineID) + "_" + strconv.Itoa(j)
				if err := writer.Write([]uint16{uint16(routineID), uint16(j)}); err != nil {
					t.Error(err)
				}
			}
		}(writer, i)
	}

	wg.Wait()
}
