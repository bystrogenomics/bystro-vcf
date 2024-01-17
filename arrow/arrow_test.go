package arrow

import (
	"fmt"
	"log"
	"os"
	"reflect"
	"sort"
	"sync"
	"testing"

	"github.com/apache/arrow/go/v14/arrow"
	"github.com/apache/arrow/go/v14/arrow/array"
	"github.com/apache/arrow/go/v14/arrow/ipc"
)

func sortRows(rows [][]interface{}) {
	sort.Slice(rows, func(i, j int) bool {
		for col := 0; col < len(rows[i]); col++ {
			val1 := rows[i][col]
			val2 := rows[j][col]
			if val1 == nil || val2 == nil {
				if (val1 == nil && val2 == nil) || val1 != nil {
					return false
				}

				return true
			}

			switch val1 := rows[i][col].(type) {
			case int:
				val2 := rows[j][col].(int)
				if val1 != val2 {
					return val1 < val2
				}
			case int8:
				val2 := rows[j][col].(int8)
				if val1 != val2 {
					return val1 < val2
				}
			case int16:
				val2 := rows[j][col].(int16)
				if val1 != val2 {
					return val1 < val2
				}
			case int32:
				val2 := rows[j][col].(int32)
				if val1 != val2 {
					return val1 < val2
				}
			case int64:
				val2 := rows[j][col].(int64)
				if val1 != val2 {
					return val1 < val2
				}
			case uint:
				val2 := rows[j][col].(uint)
				if val1 != val2 {
					return val1 < val2
				}
			case uint8:
				val2 := rows[j][col].(uint8)
				if val1 != val2 {
					return val1 < val2
				}
			case uint16:
				val2 := rows[j][col].(uint16)
				if val1 != val2 {
					return val1 < val2
				}
			case uint32:
				val2 := rows[j][col].(uint32)
				if val1 != val2 {
					return val1 < val2
				}
			case uint64:
				val2 := rows[j][col].(uint64)
				if val1 != val2 {
					return val1 < val2
				}
			case float32:
				val2 := rows[j][col].(float32)
				if val1 != val2 {
					return val1 < val2
				}
			case float64:
				val2 := rows[j][col].(float64)
				if val1 != val2 {
					return val1 < val2
				}
			case string:
				val2, _ := rows[j][col].(string)
				if val1 != val2 {
					return val1 < val2
				}
			case bool:
				val2, _ := rows[j][col].(bool)
				if val1 != val2 {
					return !val1 && val2
				}
			}
		}
		return false
	})
}

func readArrowRows(filePath string) ([][]interface{}, error) {
	// Open the file for reading
	file, err := os.Open(filePath)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	// Create a new IPC reader
	reader, err := ipc.NewFileReader(file)
	if err != nil {
		return nil, err
	}
	defer reader.Close()

	var readRows [][]interface{}

	for i := 0; i < reader.NumRecords(); i++ {
		record, err := reader.Record(i)
		if err != nil {
			return nil, err
		}

		for rowIdx := 0; rowIdx < int(record.NumRows()); rowIdx++ {
			var row []interface{}
			for _, col := range record.Columns() {
				if col.Len() <= rowIdx {
					return nil, fmt.Errorf("column length (%d) is less than row index (%d)", col.Len(), rowIdx)
				}

				if col.IsNull(rowIdx) {
					row = append(row, nil)
					continue
				}

				var foundValue interface{}
				switch v := col.(type) {
				case *array.Uint8:
					foundValue = v.Value(rowIdx)
				case *array.Uint16:
					foundValue = v.Value(rowIdx)
				case *array.Uint32:
					foundValue = v.Value(rowIdx)
				case *array.Uint64:
					foundValue = v.Value(rowIdx)
				case *array.Int8:
					foundValue = v.Value(rowIdx)
				case *array.Int16:
					foundValue = v.Value(rowIdx)
				case *array.Int32:
					foundValue = v.Value(rowIdx)
				case *array.Int64:
					foundValue = v.Value(rowIdx)
				case *array.Float16:
					foundValue = v.Value(rowIdx)
				case *array.Float32:
					foundValue = v.Value(rowIdx)
				case *array.Float64:
					foundValue = v.Value(rowIdx)
				case *array.String:
					foundValue = v.Value(rowIdx)
				case *array.Boolean:
					foundValue = v.Value(rowIdx)
				default:
					return nil, fmt.Errorf("unsupported data type: %s", col.DataType())
				}

				row = append(row, foundValue)
			}

			readRows = append(readRows, row)
		}
	}

	return readRows, nil
}

func readAndVerifyArrowFile(filePath string, expectedRows [][]interface{}, sort bool) {
	readRows, err := readArrowRows(filePath)

	if err != nil {
		log.Fatal(err)
	}

	if sort {
		sortRows(readRows)
		sortRows(expectedRows)
	}

	for i, row := range readRows {
		if !reflect.DeepEqual(row, expectedRows[i]) {
			log.Fatalf("Mismatch at row %d: got %v, want %v", i, row, expectedRows[i])
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

	if err := builder.Release(); err != nil {
		t.Fatal(err)
	}

	if err := writer.Close(); err != nil {
		t.Fatal(err)
	}

	readAndVerifyArrowFile(filePath, rows, false)
}

func TestArrowWriterHandlesNullValues(t *testing.T) {
	filePath := "null_values.feather"
	fieldTypes := []arrow.DataType{arrow.PrimitiveTypes.Uint8, arrow.PrimitiveTypes.Uint16, arrow.PrimitiveTypes.Uint32, arrow.PrimitiveTypes.Uint64,
		arrow.PrimitiveTypes.Int8, arrow.PrimitiveTypes.Int16, arrow.PrimitiveTypes.Int32, arrow.PrimitiveTypes.Int64,
		arrow.PrimitiveTypes.Float32, arrow.PrimitiveTypes.Float64, arrow.BinaryTypes.String, arrow.FixedWidthTypes.Boolean}
	batchSize := 100

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

	if err := builder.Release(); err != nil {
		t.Fatal(err)
	}

	if err := writer.Close(); err != nil {
		t.Fatal(err)
	}

	readAndVerifyArrowFile(filePath, rows, false)
}

func TestArrowWriterConcurrency(t *testing.T) {
	filePath := "concurrent_output.feather"
	fieldNames := []string{"Field1", "Field2"}
	fieldTypes := []arrow.DataType{arrow.PrimitiveTypes.Uint16, arrow.PrimitiveTypes.Uint16}

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
				rowToWrite := rows[routineID*numWritesPerRoutine+j]

				if err := builder.WriteRow(rowToWrite); err != nil {
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

	readAndVerifyArrowFile(filePath, rows, true)
}
