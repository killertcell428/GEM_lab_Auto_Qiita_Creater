'use client';

interface LoadingSpinnerProps {
  size?: 'sm' | 'md' | 'lg';
  className?: string;
}

export default function LoadingSpinner({ size = 'md', className = '' }: LoadingSpinnerProps) {
  const sizeClasses = {
    sm: 'h-4 w-4 border-b',
    md: 'h-8 w-8 border-b-2',
    lg: 'h-12 w-12 border-b-4'
  };

  return (
    <div className={`flex justify-center items-center py-8 ${className}`}>
      <div className={`animate-spin rounded-full ${sizeClasses[size]} border-blue-600`} role="status" aria-label="読み込み中">
        <span className="sr-only">読み込み中...</span>
      </div>
    </div>
  );
}

