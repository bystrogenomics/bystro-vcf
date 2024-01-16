package arrow

import (
	"log"
	"os"
	"reflect"
	"sync"
	"testing"

	"github.com/apache/arrow/go/v14/arrow"
	"github.com/apache/arrow/go/v14/arrow/array"
	"github.com/apache/arrow/go/v14/arrow/ipc"
)

func readAndVerifyArrowFile(filePath string, batchSize int, rows [][]interface{}) {
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

		for colIdx, col := range record.Columns() {
			for rowIdx := 0; rowIdx < batchSize; rowIdx++ {
				originalIdx := batchSize*i + rowIdx
				if originalIdx >= len(rows) {
					break
				}

				expectedValue := rows[originalIdx][colIdx]
				var foundValue interface{}

				switch v := col.(type) {
				case *array.Uint16:
					if v.Len() <= rowIdx {
						log.Fatal("Uint16 column length does not match number of rows")
					}
					foundValue = v.Value(rowIdx)
				case *array.Float64:
					if v.Len() <= rowIdx {
						log.Fatal("Float64 column length does not match number of rows")
					}
					foundValue = v.Value(rowIdx)
				// Add cases for other types as needed
				default:
					log.Fatalf("Unsupported column type: %v", reflect.TypeOf(col))
				}

				if !reflect.DeepEqual(foundValue, expectedValue) {
					log.Fatalf("Value in column does not match value in row: found %v, expected %v", foundValue, expectedValue)
				}
			}
		}
	}
}

func TestArrowWriteRead(t *testing.T) {
	batchSize := 5
	filePath := "test_matrix.feather"

	// Define the data types for the fields
	fieldTypes := []arrow.DataType{arrow.PrimitiveTypes.Uint16, arrow.PrimitiveTypes.Uint16, arrow.PrimitiveTypes.Uint16}
	fieldNames := []string{"Field1", "Field2", "Field3"}
	rows := make([][]interface{}, 100)
	for i := range rows {
		rows[i] = []interface{}{uint16(i), uint16(i + 1), uint16(i + 2)}
	}

	writer, err := NewArrowWriter(filePath, fieldNames, fieldTypes, batchSize)
	if err != nil {
		t.Fatal(err)
	}

	for _, row := range rows {
		if err := writer.Write(row); err != nil {
			t.Fatal(err)
		}
	}

	if err := writer.Close(); err != nil {
		t.Fatal(err)
	}

	readAndVerifyArrowFile(filePath, batchSize, rows)
}

func TestArrowWriterConcurrency(t *testing.T) {
	fieldNames := []string{"Field1", "Field2"}
	fieldTypes := []arrow.DataType{arrow.PrimitiveTypes.Uint16, arrow.PrimitiveTypes.Uint16}
	writer, err := NewArrowWriter("concurrent_output.feather", fieldNames, fieldTypes, 10)
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
				row := []interface{}{uint16(routineID), uint16(j)}
				if err := writer.Write(row); err != nil {
					t.Error(err)
				}
			}
		}(writer, i)
	}

	wg.Wait()
}
