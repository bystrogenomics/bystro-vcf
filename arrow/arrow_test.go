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

// sortRows sorts a slice of slices of any.
func sortRows(rows [][]any) {
	sort.Slice(rows, func(i, j int) bool {
		for col := 0; col < len(rows[i]) && col < len(rows[j]); col++ {
			if compareInterface(rows[i][col], rows[j][col]) {
				return true
			} else if compareInterface(rows[j][col], rows[i][col]) {
				return false
			}
		}
		return false
	})
}

// compareInterface compares two any values and returns true if the first is less than the second.
func compareInterface(a, b any) bool {
	// Handle nil values
	if a == nil || b == nil {
		return a == nil && b != nil
	}

	// Use reflection to compare values
	av := reflect.ValueOf(a)
	bv := reflect.ValueOf(b)

	// Ensure types are comparable
	if av.Type() != bv.Type() {
		panic("mismatched types")
	}

	switch av.Kind() {
	case reflect.Int, reflect.Int8, reflect.Int16, reflect.Int32, reflect.Int64:
		return av.Int() < bv.Int()
	case reflect.Uint, reflect.Uint8, reflect.Uint16, reflect.Uint32, reflect.Uint64:
		return av.Uint() < bv.Uint()
	case reflect.Float32, reflect.Float64:
		return av.Float() < bv.Float()
	case reflect.String:
		return av.String() < bv.String()
	case reflect.Bool:
		return !av.Bool() && bv.Bool()
	default:
		panic("unsupported type")
	}
}

func readArrowRows(filePath string) ([][]any, error) {
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

	var readRows [][]any

	for i := 0; i < reader.NumRecords(); i++ {
		record, err := reader.Record(i)
		if err != nil {
			return nil, err
		}

		for rowIdx := 0; rowIdx < int(record.NumRows()); rowIdx++ {
			var row []any
			for _, col := range record.Columns() {
				if col.Len() <= rowIdx {
					return nil, fmt.Errorf("column length (%d) is less than row index (%d)", col.Len(), rowIdx)
				}

				if col.IsNull(rowIdx) {
					row = append(row, nil)
					continue
				}

				var foundValue any
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

func readAndVerifyArrowFile(filePath string, expectedRows [][]any, sort bool) {
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
	rows := make([][]any, 100)
	for i := range rows {
		rows[i] = []any{uint16(i), uint16(i + 1), uint16(i + 2)}
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

	rows := make([][]any, len(fieldTypes))
	fieldNames := make([]string, len(fieldTypes))
	for i := range rows {
		rows[i] = []any{nil, nil, nil, nil, nil, nil, nil, nil, nil, nil, nil, nil}
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

	rows := make([][]any, numGoroutines*numWritesPerRoutine)
	for i := 0; i < numGoroutines; i++ {
		for j := 0; j < numWritesPerRoutine; j++ {
			rows[i*numWritesPerRoutine+j] = []any{uint16(i), uint16(j)}
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
