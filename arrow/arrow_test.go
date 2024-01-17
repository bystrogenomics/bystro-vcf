package arrow

import (
	"fmt"
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

	writer, err := NewArrowWriter(filePath, fieldNames, fieldTypes)
	if err != nil {
		t.Fatal(err)
	}

	builder, err := NewArrowRowBuilder(writer, batchSize)
	if err != nil {
		t.Fatal(err)
	}

	for _, row := range rows {
		if err := builder.WriteRow(row); err != nil {
			t.Fatal(err)
		}
	}

	if err := writer.Close(); err != nil {
		t.Fatal(err)
	}

	if err := builder.Release(); err != nil {
		t.Fatal(err)
	}

	readAndVerifyArrowFile(filePath, batchSize, rows)
}

func TestArrowWriterHandlesNullValues(t *testing.T) {
	filePath := "null_values.feather"
	fieldTypes := []arrow.DataType{arrow.PrimitiveTypes.Uint8, arrow.PrimitiveTypes.Uint16, arrow.PrimitiveTypes.Uint32, arrow.PrimitiveTypes.Uint64,
		arrow.PrimitiveTypes.Int8, arrow.PrimitiveTypes.Int16, arrow.PrimitiveTypes.Int32, arrow.PrimitiveTypes.Int64,
		arrow.PrimitiveTypes.Float32, arrow.PrimitiveTypes.Float64, arrow.BinaryTypes.String, arrow.FixedWidthTypes.Boolean}
	batchSize := 10

	rows := make([][]interface{}, len(fieldTypes))
	fieldNames := make([]string, len(fieldTypes))
	for i := range rows {
		rows[i] = []interface{}{nil, nil, nil, nil, nil, nil, nil, nil, nil, nil, nil, nil}

		switch fieldTypes[i] {
		case arrow.PrimitiveTypes.Uint8:
			rows[i][i] = uint8(i)
		case arrow.PrimitiveTypes.Uint16:
			rows[i][i] = uint16(i)
		case arrow.PrimitiveTypes.Uint32:
			rows[i][i] = uint32(i)
		case arrow.PrimitiveTypes.Uint64:
			rows[i][i] = uint64(i)
		case arrow.PrimitiveTypes.Int8:
			rows[i][i] = int8(i)
		case arrow.PrimitiveTypes.Int16:
			rows[i][i] = int16(i)
		case arrow.PrimitiveTypes.Int32:
			rows[i][i] = int32(i)
		case arrow.PrimitiveTypes.Int64:
			rows[i][i] = int64(i)
		case arrow.PrimitiveTypes.Float32:
			rows[i][i] = float32(i)
		case arrow.PrimitiveTypes.Float64:
			rows[i][i] = float64(i)
		case arrow.BinaryTypes.String:
			rows[i][i] = fmt.Sprintf("Field %d", i)
		case arrow.FixedWidthTypes.Boolean:
			rows[i][i] = true
		}
		fieldNames[i] = fmt.Sprintf("Field %d", i)
	}

	writer, err := NewArrowWriter(filePath, fieldNames, fieldTypes)
	if err != nil {
		t.Fatal(err)
	}

	builder, err := NewArrowRowBuilder(writer, batchSize)
	if err != nil {
		t.Fatal(err)
	}

	for _, row := range rows {
		if err := builder.WriteRow(row); err != nil {
			t.Fatal(err)
		}
	}

	if err := writer.Close(); err != nil {
		t.Fatal(err)
	}

	if err := builder.Release(); err != nil {
		t.Fatal(err)
	}

	readAndVerifyArrowFile(filePath, batchSize, rows)
}

func TestArrowWriterConcurrency(t *testing.T) {
	filePath := "concurrent_output.feather"
	fieldNames := []string{"Field1", "Field2"}
	fieldTypes := []arrow.DataType{arrow.PrimitiveTypes.Uint16, arrow.PrimitiveTypes.Uint16}
	batchSize := 10

	writer, err := NewArrowWriter(filePath, fieldNames, fieldTypes)
	if err != nil {
		t.Fatal(err)
	}

	var wg sync.WaitGroup

	numGoroutines := 5
	numWritesPerRoutine := 10

	rows := make([][]interface{}, numGoroutines*numWritesPerRoutine)
	for i := 0; i < numGoroutines; i++ {
		for j := 0; j < numWritesPerRoutine; j++ {
			rows[i*numWritesPerRoutine+j] = []interface{}{uint16(i), uint16(j)}
		}
	}

	for i := 0; i < numGoroutines; i++ {
		wg.Add(1)
		go func(writer *ArrowWriter, routineID int) {
			batchSize := 1 + routineID%10
			builder, err := NewArrowRowBuilder(writer, batchSize)
			if err != nil {
				t.Fatal(err)
			}

			defer wg.Done()
			for j := 0; j < numWritesPerRoutine; j++ {
				if err := builder.WriteRow(rows[routineID*numWritesPerRoutine+j]); err != nil {
					t.Fatal(err)
				}
			}

			if err := builder.Release(); err != nil {
				t.Fatal(err)
			}
		}(writer, i)
	}

	wg.Wait()

	if err := writer.Close(); err != nil {
		t.Fatal(err)
	}

	readAndVerifyArrowFile(filePath, batchSize, rows)
}
