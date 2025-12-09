'use client';

interface ErrorMessageProps {
  message: string;
  onDismiss?: () => void;
  className?: string;
}

export default function ErrorMessage({ message, onDismiss, className = '' }: ErrorMessageProps) {
  return (
    <div className={`p-4 bg-red-100 border border-red-400 text-red-700 rounded ${className}`}>
      <div className="flex justify-between items-center">
        <span>{message}</span>
        {onDismiss && (
          <button
            onClick={onDismiss}
            className="ml-4 text-red-700 hover:text-red-900 font-bold"
            aria-label="エラーメッセージを閉じる"
          >
            ×
          </button>
        )}
      </div>
    </div>
  );
}

