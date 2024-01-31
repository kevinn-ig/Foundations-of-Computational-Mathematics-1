        # Plot Newton's Method
        plt.figure(figsize=(8, 6))
        plt.plot(range(1, len(error_list)+1), error_list, marker='o', linestyle='-', label="Error")
        plt.xlabel('Iteration Number')
        plt.ylabel('Error')
        plt.title("Newton's Method - Convergence Rate")
        plt.legend()
        plt.grid(True)
        plt.show()
        output_directory = '/Users/kevin_smith/Desktop/FSU_Relevant_Stuff/fall_2023/FCM1/homework/programming/hw4/figures/prob2'
        plt.savefig(os.path.join(output_directory, 'quad_conv.png'))
